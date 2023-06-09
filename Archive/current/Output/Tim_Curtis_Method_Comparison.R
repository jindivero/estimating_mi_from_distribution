library(dplyr)
library(tidyr)
library(sdmTMB)

##Set ggplot themes
theme_set(theme_bw(base_size = 15))
theme_update(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#load simulated data
data_sims <- readRDS("~/Dropbox/Mac/Documents/GitHub/thresholds_mi_distribution/thresholds_mi_distribution/data_sims.rds")

#Add calculation of MI with 1.2
#function
add_mi <- function(x){
  mi_1.2<- x[,"po2"]*exp(0.3* x[,"invtemp"])
  mi_1.2_sc <- scale(mi_1.2)
  x <- cbind(x, mi_1.2, mi_1.2_sc)
  return(x)
}

#Apply
data_sims_d1_test <- lapply(data_sims_d1, add_mi)
data_sims_CV_d1_test <- lapply(data_sims_CV_d1, add_mi)

##Fit model
run_tim  <- function(x, mesh){
  m2 <- try(sdmTMB(observed ~ 1+as.factor(year)+breakpt(mi_0.3_sc)+log_depth_sc+log_depth_sc2, 
                   data = x, 
                   spatial = "on",
                   spatiotemporal="off",
                   mesh=mesh[[2]],
                   family =tweedie(link="log"),
                   control = sdmTMBcontrol(
                     newton_loops = 2)))
  try(tidy(m2))
  try(return(m2))
}


run_tim2  <- function(x, mesh){
m2 <- try(sdmTMB(observed ~ 1+as.factor(year)+breakpt(mi_1.2_sc)+log_depth_sc+log_depth_sc2, 
                 data = x, 
                 spatial = "on",
                 spatiotemporal="off",
                 mesh=mesh[[2]],
                 family =tweedie(link="log"),
                 control = sdmTMBcontrol(
                   newton_loops = 2)))
try(tidy(m2))
try(return(m2))
}

fits_tim <- mapply(run_tim, data_sims_d1_test, data_sims, SIMPLIFY=F)
fits_tim2 <- mapply(run_tim2, data_sims_d1_test, data_sims, SIMPLIFY=F)

#Extract pars
extract_pars <- function(x){
  if(!is.character(x)){
    par_estimates <- as.data.frame(tidy(x, conf.int = TRUE, effects="fixed"))
    par_estimates_rand <- as.data.frame(tidy(x, conf.int = TRUE, effects="ran_pars"))
    par_estimates <- bind_rows(par_estimates, par_estimates_rand)
    return(par_estimates)
  }
  if(is.character(x)){
    return(NA)
  }
}

pars_tim1 <- lapply(fits_tim, extract_pars)
pars_tim2 <- lapply(fits_tim2, extract_pars)

clean_pars <- function(pars, fits){
  names(pars) <- c(1:length(fits))
  #Remove models with errors
  pars <- keep(pars, function(x) !is.logical(x))
  #Combine into single dataframe, with column of simulation number
  pars <- bind_rows(pars,.id="id")
  return(pars)
}

##Apply to each
pars_tim1 <- clean_pars(pars_tim1, fits=fits_tim)
pars_tim2 <- clean_pars(pars_tim2, fits=fits_tim2)

##Calculate RMSE
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
  pars <- pars[c(1:3)]
  pars <- pivot_wider(pars, id_cols=c(id), names_from=term, values_from=estimate)
}
#Make wide
pars_wide_tim1 <- make_pars_wide(pars_tim1)
pars_wide_tim2 <- make_pars_wide(pars_tim2)

#Extract tweedie parameters (phi's and p's)
phis1 <- as.matrix(pars_wide_tim1[, 13])
phis2<- as.matrix(pars_wide_tim2[, 13])
ps1<- as.matrix(pars_wide_tim1[, 15])
ps2 <- as.matrix(pars_wide_tim2[, 15])

###Make predictions from model parameters
##For on original simulated data
#Remove mesh from data list
data_sims_d <- flatten(data_sims)
data_sims_d1 <- keep(.x=data_sims_d, .p=is.data.frame)

#Apply to all simulations
preds_tim1<- mapply(FUN=predict_sims, fits_tim, data_sims_CV_d1, ps1, phis1, SIMPLIFY=F)
preds_tim2<- mapply(FUN=predict_sims, fits_tim2, data_sims_CV_d1_test, ps2, phis2, SIMPLIFY=F)

###Calculate root mean squared error
##Create function to calculate RMSE 
calculate_rmse <- function(data_observed, data_predicted, column_obs, column_preds){
  observed <- as.numeric(data_observed[ , column_obs])
  predicted <- as.numeric(data_predicted[, column_preds])
  rmse1 <- rmse(actual=observed, predicted=predicted)
  return(rmse1)
}

##Apply for cross-validation
rmse_tim1<- mapply(FUN=calculate_rmse, data_sims_CV_d1, preds_tim1, "observed","pred2", SIMPLIFY=F)
rmse_tim2<- mapply(FUN=calculate_rmse, data_sims_CV_d1_test, preds_tim2, "observed","pred2", SIMPLIFY=F)

##Combine and compare
##Combine into one
rmse_all <- as.data.frame(unlist(rmse_CV1))
rmse_all$model2 <- as.vector(unlist(rmse_CV2))
rmse_all$model3 <- as.vector(unlist(rmse_CV3))
rmse_all$tim1  <- as.vector(unlist(rmse_tim1))
rmse_all$tim2 <- as.vector(unlist(rmse_tim2))

colnames(rmse_all) <- c("Data Limited", "Data Rich", "Unusual Case", "Tim 0.3", "Tim 1.2")
rmse_all$sim <- 1:100
#Pivot long
rmse_all <- pivot_longer(rmse_all, 1:5, names_to="Model")

ggplot(rmse_all, aes(x=Model, y=value, group=Model, color=Model))+geom_boxplot()+theme(legend.position="none")


#Curtis method
####Compare to presence-absence (Curtis method)
#Calculate 
##Calculate area under receiving operating curve for presence/absence
calculate_mi_quant_d<- function(dat,quant, type){
  if (type==1) {
    dat$pres <- ifelse(dat$y_0.3>0,1,0)
    test_dat <- arrange(dat, mi_0.3)
    test_dat$csum <- cumsum(test_dat$pres)
    test_dat$csum_prop <- test_dat$csum/sum(test_dat$pres)
    crit <- test_dat[which.min(abs(quant-test_dat$csum_prop)),]
    crit <- crit$mi_0.3
  }
  if (type==2) {
    dat$pres <- ifelse(dat$pred2>0,1,0)
    test_dat <- arrange(dat, mi2)
    test_dat$csum <- cumsum(test_dat$pres)
    test_dat$csum_prop <- test_dat$csum/sum(test_dat$pres)
    crit <- test_dat[which.min(abs(quant-test_dat$csum_prop)),]$mi2
  }
  return(crit)
}

dat_crits_d <- lapply(data_sims_d, calculate_mi_quant_d,quant=0.05, type=1)


#Create a presence-absence column
auc(response, predicted)



#