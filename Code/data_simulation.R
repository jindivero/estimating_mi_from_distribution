#Packages
install.packages("remotes")
library(remotes)
install.packages("devtools")
library(devtools)
install.packages("pkgbuild")
library(pkgbuild)
remotes::install_github("pbs-assess/sdmTMB", dependencies = TRUE, force=TRUE, ref="mi")
library(sdmTMB)
library(here)
library(mvtnorm)
library(mgcv)
library(dplyr)
library(devtools)
library(zoo)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(MASS)
install.packages("ggpubr")
library(ggpubr)
library(scales)
library(visreg)
library(ggeffects)
library(stringr)
library(TMB)

##Set ggplot themes
theme_set(theme_bw(base_size = 15))
theme_update(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

####Simulating data from sablefish trawl data
###Data
here("thresholds_mi_distribution")
setwd("~/Dropbox/Mac/Documents/GitHub/thresholds_mi_distribution/thresholds_mi_distribution")
dat <- readRDS("data/data_sablefish2.rds")

###Environmental data simulation
simulate_environmental <- function(ndata, dat, po2_im, temp_lim, depth_lim){
  ##Set seed
  seed <- sample(1:1000, 1)
  set.seed(seed)
  ##Simulate oxygen and temperature data
  #Get variance covariances of depth, temperature and po2
  env.dat <- dplyr::select(dat, "po2", "depth", "temp")
  var.covar <- var(log(env.dat))
  env.bar <- colMeans(log(env.dat))
  #Simulate from MVN
  log.env <- matrix(nrow=ndata, ncol=3, NA)
  colnames(log.env) <- c("po2", "depth", "temp")
  #Random generation of environmental data, with constraints
  for (i in 1:ndata){
    log.env[i,] <- rmvnorm(n = 1, mean=(env.bar - diag(var.covar) / 2), sigma=(var.covar), method="eigen")
    while(log.env[i,1]>log(po2_lim) | log.env[i,2]>log(depth_lim)| log.env[i,3]>log(temp_lim)){
      log.env[i,] <- rmvnorm(n = 1, mean=(env.bar - diag(var.covar) / 2), sigma=(var.covar), method="eigen")
    }
  }
  #Exponentiate
  env <- as.data.frame(exp(log.env))
  ##Metabolic Index
  #Constants
  #Ao <- 1.81
  #n <- -0.12
  #B = 1000 # size in grams, roughly average
  #avgbn <-B^n  # this works for the adult class (0.5 - 6 kg).  for the large adult, adjust
  kelvin = 273.15
  boltz = 0.000086173324
  tref <- 12
  #Calculate inverse temp
  env$invtemp <- (1 / boltz)  * ( 1 / (env$temp + 273.15) - 1 / (tref + 273.15))
  #Calculate MI
  env$mi_0.3 = env$po2*exp(0.3* env$invtemp)
  
  ##Add variables to dataset
  #La/lon
  env$lat <- dat$latitude
  env$lon <- dat$longitude
  #Year
  env$year <- dat$year
  ##Log and scale depth
  env$log_depth <- log(env$depth)
  env$log_depth_sc <- scale(env$log_depth)
  env$log_depth_sc2 <- env$log_depth_sc ^ 2
  ##Scale min and oxygen, just in casee
  env$mi_0.3_sc <- scale(env$mi_0.3)
  env$po2_sc <- scale(env$po2)
  #Seed
  env$seed <- seed
  return(env)
}

##Run environmental data simulation
n <- 100 #Number of simulations
ndata <- nrow(dat) #Number of data points
#Environmental data limits
po2_lim <- 35
temp_lim <-15
depth_lim <- 1500
#Simulate data
env_sims <- map(seq_len(n), ~simulate_environmental(ndata, dat, po2_lim, temp_lim,depth_lim)) 

##Fish density distribution simulation
simulate_fish <- function(env){
  ##Simulate fish density
  # set up parameters
  x50 <- 2
  delta <- 2 # how much bigger is MI at 95% of logistic compared to 50% of logistic
  b_years <- c(0.1,0,0,0.2,0.2)
  beta1 <- 1 #depth
  beta2 <- -0.25 #depth^2
  beta3 <- 4 # maximum effect of MI: means that catch above threshold is at most exp(beta3) times larger than at MI = 0
  beta0 <- 5 # log mean catch when at average depth and MI = exceeds threshold
  phi <- 7 #estimated at 20 in real sablefish data
  p <- 1.5
  range <- 0.3
  sigma_O <- 0.5
  Eo <- 0.3
  
  env$lon <- env$lon/1000
  env$lat <- env$lat/1000
  mesh <- make_mesh(env, xy_cols = c("lon", "lat"), cutoff = .015)
  data <- sdmTMB_simulate(formula=~1+as.factor(year)+logistic(mi)+log_depth_sc+log_depth_sc2,
                          data=env,
                          family=tweedie(link="log"),
                          tweedie_p=p,
                          phi=phi,
                          range=range,
                          sigma_O=sigma_O,
                          sigma_E=NULL,
                          mesh=mesh,
                          threshold_coefs=c(x50, delta, beta3, Eo),
                          B=c(beta0, b_years, beta1, beta2))
  #Add other variables back to dataset
  data$mi_0.3 <- env$mi_0.3
  data$po2 <- env$po2
  data$invtemp <- env$invtemp
  data$temp <- env$temp
  data$depth <- env$depth
  data$log_depth_sc <- as.vector(env$log_depth_sc)
  data$log_depth_sc2 <- as.vector(env$log_depth_sc2)
  data$mi_0.3_sc <- env$mi_0.3_sc
  data$year <- env$year
  data$seed <- env$seed

  data_complete <- list(data, mesh)
  return(data_complete)
}

##Run data simulations
data_sims <- lapply(env_sims, simulate_fish)
##Run again to create cross-validation data
data_sims_CV <- lapply(env_sims, simulate_fish)

###Sanity checks on simulated data

##Smoothed kernel density plot
#Extract generated data and remove mesh
data_sims_d <- flatten(data_sims)
data_sims_d1 <- keep(.x=data_sims_d, .p=is.data.frame)
data_sims_d <- bind_rows(data_sims_d1)
#Convert units of data to smaller
dat1 <- as.data.frame(dat$cpue_kg_km2/10)
colnames(dat1) <- "observed"
dat1$seed <- "Observations"
data_sims_d$seed <- as.factor(data_sims_d$seed)
data_sims_d <- bind_rows(data_sims_d,dat1)
ggplot(data_sims_d, aes(observed, color=seed, group=seed))+geom_density()+theme(legend.position="none")

#Compared to just 10 simulated datasets
data_sims_subset <- data_sims_d1[1:10]
data_sims_subset <- bind_rows(data_sims_subset)
ggplot()+geom_density(data=data_sims_subset, aes(observed, color=seed, group=seed))+theme(legend.position="none")+geom_density(data=dat, aes(cpue_kg_km2), color="red")

##Compare variance-covariance matrix
#Real data
env.dat <- dplyr::select(dat, "po2", "depth", "temp")
var.covar <- var(log(env.dat))
var.covar <- as.data.frame(as.vector(var.covar))
var.covar$type <- c("o2-o2", "o2-depth", "o2-temp", "depth-o2", "depth-depth", "depth-temp", "temp-o2", "temp-depth", "temp-temp")
colnames(var.covar)[1] <- "value"

#Simulated data
v_cov <- function(data){
  env.dat <- dplyr::select(data, "po2", "depth", "temp")
  var.covar <- var(log(env.dat))
  var.covar <- as.data.frame(as.vector(var.covar))
  rownames(var.covar) <- c("o2-o2", "o2-depth", "o2-temp", "depth-o2", "depth-depth", "depth-temp", "temp-o2", "temp-depth", "temp-temp")
  var.covar <- as.data.frame(t(as.matrix(var.covar)))
  return(var.covar)
}

v_covs <- lapply(env_sims, v_cov)
v_covs <- bind_rows(v_covs)
v_covs$sim <- 1:100
v_covs <- pivot_longer(v_covs, 1:9, names_to="type")

#Plot
ggplot(v_covs, aes(y=value, x=type))+geom_boxplot()+facet_wrap("type", scales="free")+geom_hline(data = var.covar, aes(yintercept = value, color="blue"), size=1.2)
ggplot(v_covs, aes(y=value, x=type))+geom_boxplot()+facet_wrap("type", scales="free")

#Save simulated data
saveRDS(data_sims, "data_sims.rds")
saveRDS(data_sims_CV, "data_sims_CV.rds")
