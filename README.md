# Integrated Tagging Model: Worked Example
 
Companion code for:
 
> Hyman, A. C., Collins, A. B., Ramsay, C., Wilms, S., Allen, M. S., Barbieri, L., & Frazer, T. K. (2026). Integrating conventional tagging and acoustic telemetry improves estimates of post-release survival for a highly targeted reef fish. *ICES Journal of Marine Science*, 83(11), fsag218. https://doi.org/10.1093/icesjms/fsag218
 
This repository contains a step-by-step worked example of the integrated Bayesian tagging model developed in the paper. The model combines conventional tag-return data with acoustic telemetry to estimate **absolute post-release survival** as a function of capture-and-release covariates. The example uses **simulated data** to show how to implement the model, interpret its output, and check fit. It does not reproduce the published analysis of gag (*Mycteroperca microlepis*). The model framework here also does not make the assumption that spatial information on fishing effort is known, and therefore differs slightly from the final framework presented in the published version.
 
## Why integrate the two data streams?
 
Conventional tagging is cheap and can be deployed at scale across a wide range of environmental conditions. Recapture data alone, however, cannot separate post-release survival ($\psi$) from tag reporting probability ($\xi$): the likelihood depends only on their product, so a fish that died after release looks the same as one that survived but was never reported.
 
Acoustic telemetry observes post-release fate directly, independent of angler encounter and reporting. Sharing the survival coefficients ($\beta$) across both data streams anchors the survival surface with the telemetry data. That makes $\xi$ separately estimable, while the conventional tags provide the sample size and covariate coverage.
 
## Repository contents
 
| File | Description |
|---|---|
| `Integrated_Tagging_Model_Worked_Example.pdf` | Full tutorial: model framework, worked example, and diagnostics |
| `Hyman_2026_Integrated_Tagging_Model.R` | All R code from the tutorial as a standalone script |
| `Hyman_Modified_CJS.stan` | Stan model code |
| `Simulated Conventional tag data.csv` | Simulated conventional tag releases and recaptures |
| `Simulated Acoustic tag data.csv` | Simulated acoustic telemetry survival outcomes |
| `Simulated catch data.csv` | Simulated monthly recreational catch (encounter predictor) |

 
## Model overview
 
For a conventionally tagged fish $i$ released at $t_0$, recapture at time $t$ requires surviving the immediate post-release period, remaining at risk and unencountered through each preceding interval, then being encountered and reported:
 
$$
\Pr(R_i = t) = \psi_i \prod_{s=t_{0[i]}}^{t-1}\Big(\phi\,(1-\eta_{i,s,l[i]})\Big)\,\phi\,\eta_{i,t,l[i]}\,\xi
$$
 
$$
\Pr(R_i = 0) = 1 - \sum_{t=t_{0[i]}}^{T}\Pr(R_i = t)
$$
 
with the following submodels:
 
$$
\begin{aligned}
\text{logit}(\psi_i) &= \beta_0 + \textstyle\sum_{p} \beta_p x_{i,p} && \text{post-release survival (conventional)}\\
y_j &\sim \text{Bernoulli}(\psi^\star_j), \quad \text{logit}(\psi^\star_j) = \beta_0 + \textstyle\sum_{p} \beta_p x^\star_{j,p} && \text{post-release survival (acoustic)}\\
\phi &= \exp\big(-\exp(\gamma)\big) && \text{remain at risk}\\
\eta_{i,t,l} &= 1 - \exp\Big(-\exp\big(\alpha_0 + \textstyle\sum_{q} z_{t,q}\alpha_q + \tau_l\big)\,\sigma_i\Big) && \text{encounter}\\
\text{logit}(\xi) &= \lambda && \text{reporting}
\end{aligned}
$$
 
Here $\phi$ is per-step retention in the at-risk population (losses to natural mortality, tag shedding, and emigration), $\tau_l$ are spatial random effects on encounter (non-centered), and $\sigma_i$ is an optional size-selectivity scalar (set to 1 in the example). The Stan code evaluates the product terms with a cumulative log-risk vector, which reduces the likelihood cost from $O(T^2)$ to $O(T)$ per fish. See the PDF for full derivations and assumptions.
 
## Requirements
 
- R (≥ 4.1)
- [rstan](https://mc-stan.org/rstan/) (tested with 2.32) and a working C++ toolchain
- R packages:
```r
install.packages(c("dplyr", "tidyr", "tidyverse", "lubridate", "kableExtra",
                   "knitr", "ggplot2", "ggpubr", "cowplot", "bayesplot",
                   "sp", "sf", "raster", "rstan"))
```
 
The Stan file uses the pre-2.33 array syntax (`int x[N]`). It compiles under rstan 2.32 but not under CmdStan ≥ 2.33; update declarations to `array[N] int x` if using cmdstanr.
 
## Running the example
 
1. Clone the repository and set your working directory to the repository root.
2. Run `Hyman_2026_Integrated_Tagging_Model.R` from top to bottom. The script:
   - builds a 0.25° spatial grid and assigns releases to cells,
   - discretizes release and recapture dates into monthly time steps (right-censoring fish not recaptured),
   - constructs the design matrices $X$, $X^\star$, and $Z$,
   - compiles and fits the Stan model, saving the fit to `Simulation_results.rds`,
   - reports MCMC diagnostics ($\hat{R}$, ESS, traceplots),
   - runs posterior predictive checks for both the acoustic and conventional data.
3. The fitting step can take a while. Once `Simulation_results.rds` exists, you can comment out the fitting section and start from the diagnostics.
The example fits a single chain of 2,000 iterations for speed. For real analyses, run multiple chains (e.g., 4) and check convergence on all parameters, including `sigma_tau` and the random effects.
 
## Adapting to your own data
 
Required inputs:
 
- **Acoustic telemetry:** a binary post-release survival indicator (1 = survived, 0 = died) and release covariates.
- **Conventional tagging:** release date, release location (lat/lon), release covariates, a binary recapture indicator, and recapture dates for observed returns.
- **Encounter predictors:** time-indexed covariates of fishing pressure (here, standardized monthly recreational catch).
$X$ and $X^\star$ must have identical column structure so that $\beta$ is shared between data streams. Choose the spatial grid and time step so that enough releases and recaptures fall into each cell and interval; with few recaptures, coarser discretization is recommended.
 
## Key assumptions
 
- Encountered fish are removed from the at-risk population whether or not they are reported (unreported fish are assumed harvested or released without a tag).
- $\phi$ is constant in space and time.
- Reporting probability $\xi$ is constant.
- Covariate effects on post-release survival are the same for acoustically and conventionally tagged fish.
## Citation
 
If you use this code, please cite the paper above. BibTeX:
 
```bibtex
@article{Hyman2026,
  author  = {Hyman, A. Challen and Collins, Angela B. and Ramsay, Chloe and
             Wilms, Sean and Allen, Micheal S. and Barbieri, Luiz and
             Frazer, Thomas K.},
  title   = {Integrating conventional tagging and acoustic telemetry improves
             estimates of post-release survival for a highly targeted reef fish},
  journal = {ICES Journal of Marine Science},
  year    = {2026},
  volume  = {83},
  number  = {11},
  pages   = {fsag218},
  doi     = {10.1093/icesjms/fsag218}
}
```
 
## Contact
 
A. Challen Hyman, PhD
University of South Florida, College of Marine Science / FWC Fish and Wildlife Research Institute
 
