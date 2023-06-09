library(dplyr)
library(tidyr)
library(rlang)
library(purrr)
library(mgcv)
install.packages("Metrics")
library(Metrics)
library(viridis)
install.packages("cvAUC")
library(cvAUC)
library(ggplot2)

##Set ggplot themes
theme_set(theme_bw(base_size = 15))
theme_update(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#Load model fits if needed
load("~/Dropbox/Mac/Documents/GitHub/thresholds_mi_distribution/thresholds_mi_distribution/current/output/model_fits.Rdata")

#load simulated data
data_sims <- readRDS("~/Dropbox/Mac/Documents/GitHub/thresholds_mi_distribution/thresholds_mi_distribution/current/output/data_sims.rds")

#load cross validation data
data_sims_CV <- readRDS("~/Dropbox/Mac/Documents/GitHub/thresholds_mi_distribution/thresholds_mi_distribution/current/output/data_sims_CV.rds")

#load parameter estimates
pars <- readRDS("~/Dropbox/Mac/Documents/GitHub/thresholds_mi_distribution/thresholds_mi_distribution/current/output/parameter_estimates.rds")

#####~~~~~~~~~~~~Model Predictions~~~~~~~~~
###Function to predict from each model and each dataset
predict_sims <- function(x, new_data, ps, phis){
  if(!is.character(x)){
    preds <- predict(x, newdata=new_data, type="response", return_tmb_object=F)
    #Add observation error from predictions with Tweedie parameters from model fi
    preds$pred2 <- rTweedie(preds$est, p = ps, phi = phis)
    return(preds)
  }
} 

##Get lists of tweedie parameters to feed into function
#Make function to make parameter dataframe wide
make_pars_wide <- function(pars){
  pars <- pars[c(1:3,7)]
  pars <- pivot_wider(pars, id_cols=c(id, model), names_from=term, values_from=estimate)
}

#Make wide
pars_wide <- make_pars_wide(pars)
#Save
saveRDS(pars_wide, "parameter_wide.rds")

#Extract tweedie parameters (phi's and p's)
phis1 <- as.matrix(pars_wide[c(1:100), 16])
phis2<- as.matrix(pars_wide[c(101:200), 16])
phis3 <- as.matrix(pars_wide[c(201:300), 16])
ps1<- as.matrix(pars_wide[c(1:100), 18])
ps2 <- as.matrix(pars_wide[c(101:200), 18])
ps3<- as.matrix(pars_wide[c(201:300), 18])

###Make predictions from model parameters
##For on original simulated data
#Remove mesh from data list
data_sims_d <- flatten(data_sims)
data_sims_d1 <- keep(.x=data_sims_d, .p=is.data.frame)

#Apply to all simulations
preds1<- mapply(FUN=predict_sims, fits, data_sims_d1, ps1, phis1, SIMPLIFY=F)
preds2<- mapply(FUN=predict_sims, fits2, data_sims_d1, ps2, phis2, SIMPLIFY=F)
preds3<- mapply(FUN=predict_sims, fits3, data_sims_d1, ps3, phis3, SIMPLIFY=F)

##For new data (cross-validation)
##Remove mesh from data list
data_sims_CV_d <- flatten(data_sims_CV)
data_sims_CV_d1 <- keep(.x=data_sims_CV_d, .p=is.data.frame)

##Predict on CV data for each model
preds_CV1 <- mapply(FUN=predict_sims, fits, data_sims_CV_d1, ps1, phis1, SIMPLIFY=F)
preds_CV2 <- mapply(FUN=predict_sims, fits2, data_sims_CV_d1, ps2, phis2, SIMPLIFY=F)
preds_CV3 <- mapply(FUN=predict_sims, fits3, data_sims_CV_d1, ps3, phis3, SIMPLIFY=F)

## Save all predictions
save(preds1,preds2,preds3, preds_CV1, preds_CV2, preds_CV3, file="preds.RData")
