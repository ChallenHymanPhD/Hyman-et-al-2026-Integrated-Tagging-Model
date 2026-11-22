#==========================================================================#
#              Integrated Tagging Model Worked Example                     #
#              A. Challen Hyman, PhD                                       #
#                                                                          #
#  Companion code to the worked-example PDF for Hyman et al. (2026).       #
#  Estimates post-release survival from conventional tagging and           #
#  acoustic telemetry data.                                                #
#                                                                          #
#  Required files (place in your working directory):                       #
#    - Simulated Conventional tag data.csv                                 #
#    - Simulated Acoustic tag data.csv                                     #
#    - Simulated catch data.csv                                            #
#    - Detailed_Florida_State_Boundary.shp (+ sidecar files)               #
#    - Hyman_Modified_CJS.stan                                             #
#==========================================================================#


#----------------~-------------------~------------------~------------------#
#                                Libraries                                 #
#----------------~-------------------~------------------~------------------#

# Syntax packages
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(tidyverse))
suppressMessages(library(lubridate))
suppressMessages(library(kableExtra))
suppressMessages(library(knitr))

# Visualization packages
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(cowplot))
suppressMessages(library(bayesplot))

# Spatial packages
suppressMessages(library(sp))
suppressMessages(library(sf))
suppressMessages(library(raster))

# Modeling
suppressMessages(library(rstan))

#----------------~-------------------~------------------~------------------#
#                                Housekeeping                              #
#----------------~-------------------~------------------~------------------#
# Clear out old files in R
#rm(list=ls(all=TRUE))

# set seed for reproducibility
set.seed(1234)

# "Not in" function
"%nin%" <- Negate("%in%")

# Set ggplot theme for plotting
Custom_theme <- function(){
  theme_bw(base_family = 'serif')%+replace%
    theme(axis.text = element_text(size = 16, color = 1),
          strip.text.x = element_text(size = 14,
                                      margin = margin(0.25,0,0.25,0, "cm")),
          strip.text.y = element_text(size = 14,
                                      margin = margin(0, 0.25, 0, 0.25, "cm"),
                                      angle = 270),
          axis.title = element_text(size = 18),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 12))}
theme_set(Custom_theme())

#----------------~-------------------~------------------~------------------#
#                                Load data                                 #
#----------------~-------------------~------------------~------------------#

# Set working directory (set to wherever you store the files)
# setwd("path/to/your/files")
Conventional_data <- read.csv("Simulated Conventional tag data.csv")
Acoustic_data <- read.csv("Simulated Acoustic tag data.csv")
Catch_data <- read.csv("Simulated catch data.csv")


#-----------------~-------------------~------------------~-----------------#
#                               Create grid                                #
#-----------------~-------------------~------------------~-----------------#

# Define a regular lat-lon grid using expand.grid.
# This creates all combinations of longitude and latitude values.
Mygrid <- expand.grid(
  Lon = seq(-84, -83, by = 0.25), # longitude sequence (resolution = 0.25°)
  Lat = seq(27, 28, by = 0.25)+0.25)# latitude sequence shifted by +0.25°

# Convert the data.frame to a SpatialPointsDataFrame (sp package format)
coordinates(Mygrid)<-~ Lon + Lat

# Assign coordinate reference system (CRS) as WGS84 geographic coordinates
proj4string(Mygrid) <-
  CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

# Convert grid points to raster (intermediate step for polygon creation)
grid=rasterFromXYZ(Mygrid)

# Convert raster cells to polygons (each cell becomes a spatial unit)
grid=rasterToPolygons(grid)

# Re-assign CRS to ensure spatial consistency after conversion
proj4string(grid) <-
  CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

# Convert to sf object for modern spatial workflows
grid=st_as_sf(grid)

# Extract centroid coordinates of each polygon (cell center)
centroids <- st_coordinates(st_centroid(grid))

# Store centroid latitude and longitude as attributes
grid$Lat <- as.numeric(centroids[, 2])
grid$Lon <- as.numeric(centroids[, 1])

# Create a unique identifier for each grid cell
grid$layer <- seq_len(nrow(grid))


#-----------------~-------------------~------------------~-----------------#
#                             Assign grid values                           #
#-----------------~-------------------~------------------~-----------------#

# Convert conventional dataset (fish release locations) to a spatial object
# using longitude and latitude coordinates.
coordinates(Conventional_data)<-~ Lon + Lat

# Assign WGS84 geographic coordinate reference system
proj4string(Conventional_data) <-
  CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

# Convert to sf format for compatibility with simple features workflows
Conventional_data <- st_as_sf(Conventional_data)

# Spatial join: assign each fish to the grid cell it falls within
# st_within ensures points are matched to polygon containing them
Conventional_data <- st_join(Conventional_data, grid[,1], join = st_within)

# ----- Construct spatial indexing table ----- #
# Aggregate fish-level data to grid-cell level
# Here we compute mean recapture rate per grid cell
Spatial_grid <- Conventional_data %>%
  group_by(layer) %>%                 # group by grid cell ID
  summarise(Prop = mean(Recap)) %>%   # summarize response variable
  ungroup()

# Create a compact integer index for model use
Spatial_grid <- Spatial_grid %>%
  mutate(ID = seq_len(nrow(.)))

# ----- Map grid-cell ID back onto individual fish ----- #
# Assign each fish a numeric grid index
# (required for modeling frameworks like Stan)
Conventional_data$ID <-
  Spatial_grid$ID[match(Conventional_data$layer, Spatial_grid$layer)]


# ----- Visual check: plot grid and fish locations ----- #

# Florida detailed county
Florida <- st_read("Detailed_Florida_State_Boundary.shp",quiet = TRUE)

# Plot results
ggplot()+
  geom_sf(data = Florida)+
  geom_sf(data = grid, aes(fill = as.factor(layer)), alpha = 0.3)+
  geom_sf(data = Conventional_data, aes(fill = as.factor(layer)), shape = 21, col = 1)+
  coord_sf(xlim = st_bbox(Conventional_data)[c(1,3)]+c(-0.0,0.5), ylim = c(st_bbox(Conventional_data)[c(2,4)]+ c(-0.1,0)))+
  scale_fill_viridis_d()+labs(fill = "Cell ID")+theme(legend.position = 'top')+
  guides(fill = guide_legend(nrow = 2))


#-----------------~-------------------~------------------~-----------------#
#                             Temporal indexing                            #
#-----------------~-------------------~------------------~-----------------#

# ----- Construct Temporal Index ----- #
# Use temporally-indexed catch data as index basis
Time_index <- Catch_data

# Assign integer time index (1, 2, ..., T)
Time_index$Time <- 1:nrow(Time_index)

# Extract month labels (useful for plotting / diagnostics)
Time_index$Month <- month(Time_index$Date, label = T)

# ----- Assign release and recapture times ----- #

# Map each release date to nearest month
Conventional_data$StartDate <-
  round_date(as.Date(Conventional_data$Release_Date), unit = "month")

# Map each recapture date to nearest month
Conventional_data$EndDate <-
  round_date(as.Date(Conventional_data$Recapture_Date), unit = "month")

# Convert dates into integer time indices
Conventional_data$Start_Index <-
  Time_index$Time[match(Conventional_data$StartDate,
                        Time_index$Date)]

Conventional_data$End_Index <-
  Time_index$Time[match(Conventional_data$EndDate,
                        Time_index$Date)]

# ----- Censoring ----- #
# Fish not recaptured within study window are NA in End_Index
# These are treated as right-censored at final time step
Conventional_data <- Conventional_data%>%
  mutate(End_Index = ifelse(is.na(End_Index),
                            max(Time_index$Time),
                            End_Index)
  )

# ----- Construct recapture index ----- #

# Recap_index stores the time step of observed recapture
# Non-recaptured fish are assigned 0 (no event observed)
Conventional_data$Recap_index <-
  ifelse(Conventional_data$Recap==1,
         Conventional_data$End_Index, 0)


#-----------------~-------------------~------------------~-----------------#
#                               Design matrices                            #
#-----------------~-------------------~------------------~-----------------#

# ----- Covariate matrices (shared structure for beta parameters) ----- #
# Conventional tagging dataset design matrix
X <- model.matrix(~Depth, data = Conventional_data)

# Acoustic tagging dataset design matrix
# Note: must match X in structure so beta coefficients are
# shared across datasets
X_star <- model.matrix(~Depth, data = Acoustic_data)

# ----- Encounter matrix ----- #
Z <- model.matrix(~Catch, data = Time_index)


#-----------------~-------------------~------------------~-----------------#
#                               Fit Stan model                             #
#-----------------~-------------------~------------------~-----------------#
# NOTE: This section is slow. To skip refitting, comment it out and load
# the saved fit in the "MCMC Diagnostics" section below.

# Compile Stan model from external file
Model <- stan_model("Hyman_Modified_CJS.stan")

# ----- Simplifying assumption: equal selectivity ----- #
# Set all fish to equal catchability (sigma = 1) for simplicity
# This removes size-based selectivity from the likelihood
Conventional_data$Selectivity <- 1

# ----- Construct data list for Stan ----- #
# This list defines the full model input structure:
# - N_CT / N_AT: sample sizes for each data source
# - T: number of time steps
# - L: number of spatial locations
# - P / Q: number of covariates for survival and encounter processes
# - Index variables define movement through space/time
# - Design matrices encode covariate effects

stan_data <- list(
  N_CT = nrow(Conventional_data),     # number of conventional tags
  N_AT = nrow(Acoustic_data),         # number of acoustic tags
  
  T = max(Time_index$Time),           # total number of time steps
  L = nrow(Spatial_grid),             # number of spatial grid cells
  
  P = ncol(X),                        # number of survival covariates
  Q = ncol(Z),                        # number of encounter covariates
  
  # -- Temporal structure -- #
  t0  = Conventional_data$Start_Index,# release time index
  R   = Conventional_data$Recap_index,# recapture time index (0 = censored)
  
  # -- Spatial structure -- #
  loc = Conventional_data$ID,         # spatial grid cell index
  
  # -- Selectivity -- #
  sigma = Conventional_data$Selectivity,
  
  # -- Covariate design matrices -- #
  X = X,                              # conventional tagging covariates
  X_star = X_star,                    # acoustic tagging covariates
  Z = Z,                              # temporal encounter covariates
  
  # -- Acoustic survival observations -- #
  y = Acoustic_data$y_star            # Acoustic tag survival indicator
)


# ----- Run Model ----- #

# Use all available cores for parallel chains
options(mc.cores = parallel::detectCores())

# Fit model using MCMC sampling
# Note: single chain + 2000 iterations is for testing/efficiency only
fit <- sampling(Model, data = stan_data,
                chains = 1, iter = 2000,
                refresh = 10)

# Fix compiled model pointer for later reuse (optional workflow detail)
fit@stanmodel@dso <- new("cxxdso")

# Save fitted model object
saveRDS(fit, file = "Simulation_results.rds")


#-----------------~-------------------~------------------~-----------------#
#                             MCMC Diagnostics                             #
#-----------------~-------------------~------------------~-----------------#

# Load model
fit <- readRDS("Simulation_results.rds")

# Define parameters to extract
Params <- c("beta", "alpha", "gamma", "lambda")

# Extract summary statistics
Summary_fit <- summary(fit)$summary

# Extract rownames
rn <- rownames(Summary_fit)

# strip Stan indexing (everything after '[')
base_names <- sub("\\[.*\\]", "", rn)

# Isolate parameters of interest
keep <- base_names %in% Params

Parameter_summary <- Summary_fit[keep,]

# Traceplots of key parameters
traceplot(fit, c("beta", "alpha", "gamma", "lambda"))


#-----------------~-------------------~------------------~-----------------#
#                            Acoustic Tag PPCs                             #
#-----------------~-------------------~------------------~-----------------#

# ----- Generate posterior predictive survival draws ----- #
# Extract posterior draws
posterior <- rstan::extract(fit)

# Extract beta posterior distributions
betas <- posterior$beta

# Generate full posterior predictive distributions for psi*
psi_star <- plogis(X_star %*% t(as.matrix(betas)))

S <- 1000  # posterior predictive replicates
N <- length(Acoustic_data$y_star)

y_rep <- matrix(NA, nrow = N, ncol = S)

for (s in 1:S) {
  for (j in 1:N) {
    y_rep[j,s] <- rbinom(1, 1, psi_star[j,s])
  }
}

# ----- Summarize simulated survival ----- #
# Average overall predicted survival
sim_means <- apply(y_rep, 2, mean)

# Survival as a function of depth
df_ppc <- tibble(
  predicted = apply(y_rep, 1, mean),
  observed = Acoustic_data$y_star
)%>%
  mutate(pred_bin = cut(predicted, breaks = seq(0, 1, by = 0.01))) %>%
  group_by(pred_bin) %>%
  summarise(
    n = n(),
    obs_mean = mean(observed),
    pred_mean = mean(predicted)
  )#%>%filter(n > 10)


Y_rep_ones <- apply(y_rep, 2, sum)%>%quantile(.,c(0.1,0.5,0.9))
Y_rep_zeros <- (N-apply(y_rep, 2, sum))%>%quantile(.,c(0.1,0.5,0.9))
Y_rep_data <- rbind(Y_rep_ones, Y_rep_zeros)%>%as.data.frame()%>%
  mutate(Value = c(1,0))


# ----- Plot results ----- #
# Histogram of observed vs predicted average post-release survival
Hist_plot_AT <-ggplot(data.frame(predicted = sim_means)) +
  geom_histogram(aes(x = predicted),
                 fill = "#B3CDE0", alpha = 0.5,
                 col = "#000A39") +
  geom_vline(xintercept = mean(Acoustic_data$y_star),
             lwd = 1, col = "#000A39")+
  labs(y = "Frequency", x = "Proportion recaptured")

# Plot scatter plot with 1-1 line
Scatter_plot_AT <- ggplot(df_ppc, aes(x = pred_mean,
                                      y = obs_mean, weight = n)) +
  geom_point(aes(), fill = "#B3CDE0", alpha = 0.9,
             col = "#000A39", shape = 21, size = 5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "#000A39", lwd = 1) +
  labs(x = "Predicted survival rate", y = "Observed recapture rate")

# Plot bar plot of observed vs. predicted counts
Bar_plot_AT<- ggplot(Acoustic_data)+
  geom_bar(aes(y_star), fill = "#B3CDE0",
           alpha = 0.5, col = "#000A39")+
  geom_point(data = Y_rep_data, aes(Value, `50%`),
             alpha = 0.9, col = "#000A39", size = 5)+
  geom_linerange(data = Y_rep_data,
                 aes(Value, ymin = `10%`, ymax = `90%`),
                 alpha = 0.9, col = "#000A39", size = 2)+
  scale_x_continuous(breaks = c(0,1))+labs(y = "Frequency", x = "survival")
Top <- plot_grid(Hist_plot_AT, Bar_plot_AT,
                 ncol = 2, align = 'h', labels = c("A", "B"))
plot_grid(Top, Scatter_plot_AT,ncol = 1, labels = c("", "C"))


#-----------------~-------------------~------------------~-----------------#
#                       Simulated posterior recaptures                     #
#-----------------~-------------------~------------------~-----------------#

simulate_recapture <- function(phi_post, h_recap_post, psi_post, xi_post,
                               t0, loc, Tmax, Selectivity, n_iter) {
  # number of conventionally tagged fish
  N <- length(t0)
  
  # initialize output matrix; 0 = censored (default)
  recaps_sim <- matrix(0L, nrow = n_iter, ncol = N)
  
  # loop over posterior draws
  for (iter in 1:n_iter) {
    
    # loop over fish
    for (i in 1:N) {
      
      # Post-release survival (censored if died post-release)
      if (rbinom(1, 1, psi_post[iter, i]) == 0) next
      
      # encounter probability at each time step for fish i
      eta <- 1 - exp(-h_recap_post[iter, , loc[i]] * Selectivity[i])
      
      # loop over time steps from release to end of study
      for (t in t0[i]:Tmax) {
        
        # probability fish survived interval, was not encountered
        p_remain    <- phi_post[iter] * (1 - eta[t])
        
        # probability fish was lost (died, emigrated, or shed tag)
        p_loss      <- 1 - phi_post[iter]
        
        # probability fish survived, was encountered, and reported
        p_report    <- phi_post[iter] * eta[t] * xi_post[iter]
        
        # probability fish survived, encountered, but not reported
        p_nonreport <- phi_post[iter] * eta[t] * (1 - xi_post[iter])
        
        # draw one of four mutually exclusive outcomes
        event_draw <- sample(c("Remain",
                               "Recapture",
                               "Lost",
                               "Not reported"),
                             size = 1,
                             prob = c(p_remain, p_report,
                                      p_loss, p_nonreport))
        # Outcome 1: fish remains at risk; advance to next time step
        if (event_draw == "Remain")       next
        
        # Outcome 2: fish lost; censored
        if (event_draw == "Lost")         break
        
        # Outcome 3: fish encountered but not reported; censored
        if (event_draw == "Not reported") break
        
        # Outcome 4: fish recaptured; record time step
        if (event_draw == "Recapture"){recaps_sim[iter, i] <- t; break}
      }
      # if Tmax reached without recapture, default 0 retained
    }
  }
  # returns n_iter x N matrix of recapture times (0 = censored)
  return(recaps_sim)
}


# ----- Simulate recaptures ----- #
# Extract posterior draws from all parameters
psi_post <- posterior$psi          # post-release survival
phi_post <- posterior$phi          # remain at risk
xi_post <- plogis(posterior$lambda)# reporting rate
h_recap_post <- posterior$h_recap  # hazard of recapture

# Extract data to simulate recapture
Tmax <- max(Time_index$Time)         # max time steps
t0 <- Conventional_data$Start_Index  # release time
loc <- Conventional_data$ID          # spatial cell
Selectivity <- rep(1,nrow(Conventional_data))  # selectivity

# Run simulation function
recap_sim <- simulate_recapture(phi_post = phi_post,
                                h_recap_post = h_recap_post,
                                psi_post = psi_post,
                                xi_post = xi_post,
                                t0=t0,
                                loc = loc,
                                Tmax = Tmax,
                                Selectivity = Selectivity,
                                n_iter = 1000)

# Summarize simulated recaptures
sim_means <- apply(ceiling(recap_sim/Tmax), 1, mean)

# ----- Plot results ----- #
# Histogram of observed vs predicted
Hist_plot_CT <- ggplot(data.frame(
  predicted = sim_means
)) +
  geom_histogram(aes(x = predicted), fill = "#B3CDE0",
                 alpha = 0.5, col = "#000A39") +
  geom_vline(xintercept = mean(Conventional_data$Recap),
             lwd = 1, col = "#000A39")+
  labs(y = "Frequency", x = "Proportion recaptured")

y <- Conventional_data$Recap_index
yrep <- recap_sim[1:1000,]

Density_plot_CT <- ppc_dens_overlay(y = y, yrep = yrep) +
  ggplot2::labs(
    x = "Recaptures time step (0 = censored)",
    y = "Kernel density"
  )+theme(legend.position = 'top')

ECDF_plot_CT <- ppc_ecdf_overlay(y, yrep)+
  ggplot2::labs(
    #title = "Posterior predictive check of recapture times",
    x = "Recaptures time step (0 = censored)",
    y = "ECDF of recaptures"
  )+theme(legend.position = 'top')

df_ppc <- tibble(
  predicted = apply(ceiling(recap_sim/Tmax), 2, mean),
  observed = Conventional_data$Recap
) %>%
  mutate(pred_bin = cut(predicted, breaks = seq(0, 1, by = 0.005))) %>%
  group_by(pred_bin) %>%
  summarise(
    n = n(),
    obs_mean = mean(observed),
    pred_mean = mean(predicted)
  )%>%filter(n>10) # filters small bins to reduce noise

# Plot calibration curve
scatter_plot_CT <- ggplot(df_ppc,
                          aes(x = pred_mean, y = obs_mean, weight = n)) +
  geom_point(aes(), fill = "#B3CDE0", alpha = 0.9,
             col = "#000A39", shape = 21, size = 5) +
  #geom_smooth()+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "#000A39", lwd = 1) +
  labs(x = "Predicted recapture rate", y = "Observed recapture rate")


ggarrange(Hist_plot_CT, scatter_plot_CT,
          Density_plot_CT, ECDF_plot_CT,
          common.legend = T, labels = c("A", "B", "C", "D"),
          font.label = list(size = 20, face = "bold",
                            family = "serif"), align = "v")