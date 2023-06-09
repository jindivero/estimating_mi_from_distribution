###~~~~~~~~~~~~~~~~~~Fit model in sdmTMB~~~~~~~~~~~~~~~~
install.packages("remotes")
library(remotes)
install.packages("devtools")
library(devtools)
install.packages("pkgbuild")
library(pkgbuild)
remotes::install_github("pbs-assess/sdmTMB", dependencies = TRUE, ref="mi")
library(sdmTMB)

##Load simulated data if needed
readRDS("~/Dropbox/Mac/Documents/GitHub/thresholds_mi_distribution/thresholds_mi_distribution/data_sims.rds")

##Specify all parameters from data generation for reference
x50 <- 2
delta <- 2 # how much bigger is MI at 95% of logistic compared to 50% of logistic
beta1 <- 1
beta2 <- -0.25
beta3 <- 4 # maximum effect of MI: means that catch above threshold is at most exp(beta3) times larger than at MI = 0
beta0 <- 5 # log mean catch when at average depth and MI = exceeds threshold
b_years <- c(0.1,0,0,0.2,0.2)
phi <- 7
p <- 1.5
range <- 0.3
sigma_O <- 0.5
Eo <- 0.3

###Model 1: No priors
##Set starting parameters
#Correct values
start <- matrix(0, ncol = 1, nrow = 4)
start[1,1] <- x50
start[2,1] <- delta
start[3,1] <- beta3
start[4,1] <- Eo
#Just above zero
start2 <- matrix(0.08, ncol = 1, nrow = 4)
#Slightly farther
start3 <- matrix(0, ncol = 1, nrow = 4)
start3[1,1] <- x50*0.1
start3[2,1] <- delta*0.1
start3[3,1] <- beta3*0.1
start3[4,1] <- Eo*0.1

##Function to run model and return list of model outputs
run_sdmTMB_1 <- function(x, start) {
    m2 <- try(sdmTMB(observed ~ 1+as.factor(year)+logistic(mi)+log_depth_sc+log_depth_sc2, 
                     data = x[[1]], 
                     spatial = "on",
                     spatiotemporal="off",
                     mesh=x[[2]],
                     family =tweedie(link="log"),
                     control = sdmTMBcontrol(
                      start = list(b_threshold = start),
                     #lower = list(b_threshold = c(0, 0, 0, 0)), upper = list(b_threshold = c(Inf, Inf, Inf, Inf)),
                     newton_loops = 1)))
    try(tidy(m2))
    try(return(m2))
}

##Fit model to all simulated datasets
fits <- lapply(data_sims, run_sdmTMB_1, 
               start=start)
#fits2 <- lapply(data_sims, run_sdmTMB_1, 
                #start=start2)
#fits3 <- lapply(data_sims, run_sdmTMB_1, 
                #start=start3)

###Model 2: Prior (correct)
run_sdmTMB_2 <- function(x, start) {
  m2 <- try(sdmTMB(observed ~ 1+as.factor(year)+logistic(mi)+log_depth_sc+log_depth_sc2, 
                   data = x[[1]], 
                   spatial = "on",
                   spatiotemporal="off",
                   mesh=x[[2]],
                   family =tweedie(link="log"),
                   control = sdmTMBcontrol(
                     start = list(b_threshold = start),
                     newton_loops = 1),
                   priors=sdmTMBpriors(threshold = normal(c(NA, NA, NA, Eo), c(NA, NA, NA, 0.15)))))
  try(tidy(m2))
  try(return(m2))
}

fits2 <- lapply(data_sims, run_sdmTMB_2, 
               start=start)

###Model 3: Prior (incorrect)
run_sdmTMB_3 <- function(x, start) {
  m2 <- try(sdmTMB(observed ~ 1+as.factor(year)+logistic(mi)+log_depth_sc+log_depth_sc2, 
                   data = x[[1]], 
                   spatial = "on",
                   spatiotemporal="off",
                   mesh=x[[2]],
                   family =tweedie(link="log"),
                   control = sdmTMBcontrol(
                     start = list(b_threshold = start),
                     newton_loops = 1),
                   priors=sdmTMBpriors(threshold = normal(c(NA, NA, NA, 1.2), c(NA, NA, NA, 0.6)))))
  try(tidy(m2))
  try(return(m2))
}
start4 <- matrix(0, ncol = 1, nrow = 4)
start4[1,1] <- x50
start4[2,1] <- delta
start4[3,1] <- beta3
start4[4,1] <- 1.2

fits3 <- lapply(data_sims, run_sdmTMB_3, 
                 start=start4)

#Save models
#save(fits, fits1, fits2, fits3, fits_b, fits_c, file="model_fits.Rdata")
save(fits, fits2, fits3, file="model_fits.Rdata")
