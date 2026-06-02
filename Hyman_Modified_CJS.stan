// ---------- Data block ---------- //
// All inputs required for the likelihood and covariates
data {
  int<lower=1> N_CT;                                                            // number of conventionally tagged fish
  int<lower=1> N_AT;                                                            // number of acoustic telemetry fish
  int<lower=1> T;                                                               // total number of discrete time steps
  int<lower=1> L;                                                               // number of spatial locations
  int<lower=1> K;                                                               // number of post-release survival predictors
  int<lower=1> Q;                                                               // number of encounter predictors
  int<lower=1> t0[N_CT];                                                        // release time step for each fish
  int<lower=0> R[N_CT];                                                         // recapture time step (0 = never recaptured / censored)
  int<lower=1> loc[N_CT];                                                       // location index for each fish
  vector[N_CT] sigma;                                                           // recreational selectivity (optional, set to 1 if not available)
  matrix[N_CT,K] X;                                                             // conventional tagging predictors of post-release survival
  matrix[N_AT,K] X_star;                                                        // acoustic tagging predictors of post-release survival
  matrix[T,Q] Z;                                                                // predictors of encounter
  int<lower = 0> y[N_AT];                                                       // Fish survival in acoustic telemetry
}

parameters {
  real gamma;                                                                   // loss term
  vector[Q] alpha;                                                              // intercept and slopes for encounter based on harmonic temporal encounter covariates
  vector[K] beta;                                                               // intercept and slopes for post-release survival covariates
  vector[L] z_tau;                                                              // location-specific random effects (non-centered parameterization; tau = z_tau * sigma_tau)
  real<lower = 0> sigma_tau;                                                    // SD of the random-effects
  real lambda;                                                                  // logit-scale reporting rate
}

// ---------- Transformed parameters ---------- //
transformed parameters {
  matrix[T,L] h_recap;                                                          // hazard of re-encounter per step
  vector[L] tau;                                                                // spatial random effect on encounter
  vector[N_CT] psi = inv_logit(X*beta);                                         // conventional tag post-release survival probability
  vector[N_AT] psi_star = inv_logit(X_star*beta);                               // Acoustic taganchor point
  real phi = exp(-exp(gamma));                                                  // constant post-recovery survival
  real log_phi = -exp(gamma);                                                   // log version for computational efficiency
  vector[T] Z_alpha = Z*alpha;                                                  // pre-compute vector of temporal effects
  real log_xi = log(inv_logit(lambda));                                         // convert reporting rate to log-probability scale
  
  // Fill time × location arrays for survival and encounter
  for (l in 1:L){
    tau[l] =  z_tau[l]*sigma_tau;                                               // add random effect to re-encounter
    for (t in 1:T){
      h_recap[t,l] = exp(Z_alpha[t] + tau[l]);                                  // log-hazard of encountered in step t; eta = 1 - exp(-h_recap * sigma)
    }
  }
  
}

// ---------- Model block ---------- //
model {
  // ---- Priors ---- //
  beta ~ normal(0,2);                                                           // post-release survival predictors
  alpha ~ normal(0,2);                                                          // harmonic priors
  gamma ~ normal(0,2);                                                          // loss prior
  sigma_tau ~ normal(0,1);                                                      // SD of location random effect
  z_tau ~ normal(0,1);                                                          // standardized random effects
  lambda ~ normal(0,1.5);                                                       // Reporting rate prior in logit space
  
  // ---- acoustic tag post-release survival likelihood ---- // 
  y ~ bernoulli(psi_star);                                                     
  
  // ---- conventional tag recapture likelihood ---- // 
  for (i in 1:N_CT) {
    int t_0i = t0[i];       // release step for fish i
    int end  = R[i];        // recapture time (0 if censored)
    int l    = loc[i];      // location
    vector[T] eta = 1-exp(-h_recap[,l]*sigma[i]);                               // vector of all encounter probabilities for fish i

    // ---- cumulative risk terms ---- //
    // The cumulative log-risk vector allows efficient computation of the
    // probability of surviving without being encountered through any interval.
    // cum_log_risk[t] = sum_{k=1}^{t} [log_phi + log(1 - eta[k])]
    // so cum_log_risk[t-1] - cum_log_risk[t_0i] recovers the log of the product
    // term in the recapture equation between t_0i+1 and t-1, without recomputing
    // it from scratch for each candidate recapture time t (more efficient).
    vector[T] log_risk;                                                         // log(phi) + log(1 - eta) for each step
    vector[T] cum_log_risk;                                                     // cumulative sum of log_risk
                                                                                                                                                
    for (k in 1:T){log_risk[k] = log_phi + log1m(eta[k]);}                     
    
    cum_log_risk[1] = log_risk[1];
    for (k in 2:T){cum_log_risk[k] = cum_log_risk[k - 1] + log_risk[k];}

    if (R[i] == 0) {
      int T_follow = T;                                                         // follow-up window set to full study period with optional adjustment min(T, t_0i + max_follow_up)
      
      // -------- Fish i never recaptured (right-censored) -------- //
      // t == t_0i: recapture in the first interval
      real log_p_recapture_sum = log(psi[i]) + log(phi) + log(eta[t_0i]) + log_xi;
      
      // t > t_0i: recapture in subsequent full intervals
      for (t in (t_0i + 1):T_follow) {
        real log_term = log(psi[i]) + log(phi) + log1m(eta[t_0i]);
        if (t > t_0i + 1) {
          log_term += cum_log_risk[t - 1] - cum_log_risk[t_0i];
        }
        log_term += log_phi + log(eta[t]) + log_xi;
        log_p_recapture_sum = log_sum_exp(log_p_recapture_sum, log_term);
      }
      
      // add log probability of no recapture = log(1 - p(any recapture))
      target += log1m_exp(log_p_recapture_sum);

    } else {
      
      // -------- Fish i recaptured at 'end' -------- //
      if (end == t_0i) {
          // Recaptured in release interval
          // Fish survived first interval and was encountered before completing it
          real log_Recap_prob = log(psi[i]) + log(phi) + log(eta[t_0i]) + log_xi;
          target += log_Recap_prob;
      } else {
          real log_Recap_prob = log(psi[i]) + log(phi) + log1m(eta[t_0i]);
          for (t in (t_0i + 1):(end - 1)){
              log_Recap_prob += log_phi + log1m(eta[t]);
          }
          log_Recap_prob += log_phi + log(eta[end]) + log_xi;
          target += log_Recap_prob;
      }
    }
  }
}

