
####~~~~~~~Critical Metabolic Index Value~~~~~~~~~~~~
###With test prediction example
total <- sum(preds_test$pred2)
preds_test <- arrange(preds_test, mi2)
preds_test$csum <- cumsum(preds_test$pred2)
preds_test$csum_prop <- preds_test$csum/total

ggplot(preds_test, aes(x=mi2,y=csum_prop))+geom_line()+xlab("Metabolic Index")+ylab("Proportion of Cumulative Sum of Fish Density")

###Compare to example generated data
total <- sum(test_dat$observed)
test_dat <- arrange(test_dat, mi_0.3)
test_dat$csum <- cumsum(test_dat$observed)
test_dat$csum_prop <- test_dat$csum/total

ggplot(preds_test, aes(x=mi2,y=csum_prop))+geom_line(color="blue")+geom_line(data=test_dat,aes(x=mi_0.3,y=csum_prop))+xlab("Metabolic Index")+ylab("Proportion of Cumulative Sum of Fish Density")+geom_hline(yintercept=0.05)


##Identify critical value
crit <- preds_test[which.min(abs(0.05-preds_test$csum_prop)),]$mi2
crit_dat <- test_dat[which.min(abs(0.05-test_dat$csum_prop)),]$mi_0.3

#Plot critical metabolic index vs actual fish density

ggplot(preds_test,aes(x=mi2, y=pred2))+geom_point(color="blue")+geom_vline(xintercept=crit, color="blue")+geom_point(data=test_dat,aes(x=mi_0.3, y=observed))+geom_vline(xintercept=crit_dat)

#How sensitive is this critical metabolic index calculation to changes in Eo
Eo_q <- quantile(Eos, c(0,.25,.5,.75,1))
calculate_mi <-  function(Eo){
  preds_test$mi <- preds_test$po2*exp(Eo*preds_test$invtemp)
}

mi_q <- lapply(Eo_q, calculate_mi)
mi_q <-as.data.frame(mi_q)
mi_q$preds <- preds_test$pred2

#Function to calculate for each

total <- sum(mi_q$preds)
t <- arrange(mi_q, X0.)
t$csum <- cumsum(mi_q$preds)
t$csum_prop <- t$csum/total
crit <- t[which.min(abs(0.05-t$csum_prop)),]
crit0 <- crit$X0.
t <- arrange(mi_q, X25.)
t$csum <- cumsum(mi_q$preds)
t$csum_prop <- t$csum/total
crit <- t[which.min(abs(0.05-t$csum_prop)),]
crit25 <- crit$X25.
t <- arrange(mi_q, X50.)
t$csum <- cumsum(mi_q$preds)
t$csum_prop <- t$csum/total
crit <- t[which.min(abs(0.05-t$csum_prop)),]
crit50 <- crit$X50.
t <- arrange(mi_q, X75.)
t$csum <- cumsum(mi_q$preds)
t$csum_prop <- t$csum/total
crit <- t[which.min(abs(0.05-t$csum_prop)),]
crit75 <- crit$X75.
t <- arrange(mi_q, X100.)
t$csum <- cumsum(mi_q$preds)
t$csum_prop <- t$csum/total
crit <- t[which.min(abs(0.05-t$csum_prop)),]
crit1 <- crit$X100.

crit_sensitivity <- as.data.frame(as.vector(c(crit0, crit25, crit50,crit75,crit1)))
colnames(crit_sensitivity) <-  "MI_crit"

crit_sensitivity$Eo <- Eo_q

ggplot(crit_sensitivity, aes(x=Eo, y=MI_crit))+geom_point()+geom_line()

#Plot metabolic index
ggplot(preds_test, aes(x=lon, y=lat))+geom_point(aes(color=mi2), size=0.25)
ggplot(preds_test, aes(x=lon, y=lat))+geom_point(aes(color=mi2), size=2)+geom_point(data=subset(preds_test, pred2 > 0), aes(x=lon, y=lat), shape=21, size=.5, fill="pink")

ggplot(preds_test, aes(x=lat, y=depth_scaled))+geom_point(aes(color=mi2), size=5)+scale_y_reverse()+geom_point(data=subset(preds_test, pred2 > 0), aes(x=lat, y=depth_scaled), shape=21, size=1, color="pink")
ggplot(preds_test, aes(x=mi2, y=depth_scaled))+geom_point(aes(color=log(pred2)), size=.2)+scale_y_reverse()+geom_vline(xintercept=crit50)+scale_color_viridis()

#Plot po2 and invtemp, colored by fish density
#need to calculate critical po2 and critical-mi curve
preds_test$POcrit <- exp(-test_pars$"mi-Eo"* preds_test$invtemp)
preds_test$POcrit <- exp(-test_pars$"mi-Eo"* preds_test$invtemp)
preds_test$POcrit2 <- exp(-test_pars$"mi-Eo"* preds_test$invtemp)+crit50

ggplot(preds_test)+geom_point(aes(x=invtemp, y=po2, color=log(pred2)))+geom_line(data=preds_test, aes(x=invtemp,  y=POcrit), color="blue")+scale_x_reverse()

ggplot(subset(preds_test, pred2 > 0))+geom_point(aes(x=invtemp, y=po2, color=log(pred2)))+geom_line(data=preds_test, aes(x=invtemp,  y=POcrit), color="blue")+geom_line(data=preds_test, aes(x=invtemp,  y=POcrit2), color="blue")+scale_x_reverse()

###Make function for doing for all iterations

#for generated data
data_sims_d <- flatten(data_sims)
data_sims_d <- keep(.x=data_sims_d, .p=is.data.frame)
dat_crits <- lapply(data_sims_d, calculate_mi_quant,quant=0.05, type=1)

#for predicted data
pred_crits <- lapply(preds_mi, calculate_mi_quant,quant=0.05, type=2)

#combine
dat_crits <- unlist(dat_crits)
pred_crits <- unlist(pred_crits)
crits <- as.data.frame(cbind(dat_crits, pred_crits))
crits2 <- pivot_longer(crits, cols=1:2, names_to="sim")

ggplot(crits, aes(x=dat_crits, y=pred_crits))+geom_point()+xlab("MI-crit from generated. data")+ylab("MI-crit from predictions from model fit")

ggplot(crits2, aes(y=value, group=sim, color=sim))+geom_boxplot()




###Critical metabolic index from Deutsch et al. presence/absence
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
pred_crits_d <- lapply(preds_mi, calculate_mi_quant_d,quant=0.05, type=2)

dat_crits_d <- unlist(dat_crits_d)
pred_crits_d <- unlist(pred_crits_d)
crits_d <- as.data.frame(cbind(dat_crits_d, pred_crits_d))
crits2_d <- pivot_longer(crits_d, cols=1:2, names_to="sim")

ggplot(crits2_d, aes(y=value, group=sim, color=sim))+geom_boxplot()

###Critical metabolic index from predicted fit, rather than with observation error?
pred_crits_fit <- lapply(preds_mi, calculate_mi_quant,quant=0.05, type=3)
pred_crits_fit <- unlist(pred_crits_fit)
pred_crits_fit <- as.data.frame(cbind(pred_crits_fit))
crits <- as.data.frame(cbind(crits, pred_crits_fit))
crits2 <-  pivot_longer(crits, cols=1:3, names_to="sim")

#Combine Curtis method and using density
crits2_d$method <-"presense/absence"
crits2$method <- "density"
crits_complete <- bind_rows(crits2_d, crits2)

ggplot(crits_complete, aes(y=value, x=sim, group=sim))+geom_boxplot()+scale_x_discrete(name = "Method", labels = c("Data, Density", "Data,Presence/Absence", "Predictions,Density","Predictions, Presence/Absence", "Density no obs error"))

#Look at metabolic index threshold along depth gradient?
ggplot(preds_test, aes(x=log_depth_sc, y=mi2))+geom_point(aes(color=pred2))

#Ratio of MI to depth?

#Critical depth?
calculate_mi_quant<- function(dat,quant, type){
  if (type==1) {
    total <- sum(dat$y_0.3)
    test_dat <- arrange(dat, depth_scaled)
    test_dat$csum <- cumsum(test_dat$y_0.3)
    test_dat$csum_prop <- test_dat$csum/total
    crit <- test_dat[which.min(abs(quant-test_dat$csum_prop)),]
    crit <- crit$depth_scaled
  }
  if (type==2) {
    total <- sum(dat$pred2)
    test_dat <- arrange(dat, mi2)
    test_dat$csum <- cumsum(test_dat$pred2)
    test_dat$csum_prop <- test_dat$csum/total
    crit <- test_dat[which.min(abs(quant-test_dat$csum_prop)),]$mi2
  }
  if (type==3) {
    total <- sum(dat$pred)
    test_dat <- arrange(dat, mi2)
    test_dat$csum <- cumsum(test_dat$pred)
    test_dat$csum_prop <- test_dat$csum/total
    crit <- test_dat[which.min(abs(quant-test_dat$csum_prop)),]$mi2
  }
  return(crit)
}

#Calculating critical metabolic index
test_dat0$po2_crit <- 0.9/(exp(test_pars$"mi-Eo"* test_dat0$invtemp))

ggplot(test_dat0, aes(temp, po2, z=y_0.3)) +
  stat_summary_2d(fun=sum) +
  geom_line(data=test_dat0, aes(x=temp, y=po2_crit), color="red")









###Example of just one
##Prediction
test2 <- preds[[20]]
test <- test2[[1]]

##Predict from simulated data and model fit
test_fit <- fits[[20]]
##Simulated data
test_dat <- data_sims[[20]]
test_dat <- as.data.frame(test_dat[1])
##Estimated parameters
test_pars <- pars_wide[20,]
##Calculate metabolic index
test_dat$mi_2 <- test_dat$po2*exp(test_pars$"mi-Eo"* test_dat$invtemp)

##Prediction from simulated data, not real data
pred_dat2 <- test_dat[,c('po2', 'invtemp', 'log_depth_sc', 'log_depth_sc2', 'lat', 'lon', 'year')]
preds_test <- predict(test_fit, newdata=pred_dat2, type="response", return_tmb_object=T)
preds_test <- preds_test[[1]]
preds_test$pred <- exp(preds_test$est)
##Calculate metabolic index
preds_test$mi2 <- preds_test$po2*exp(test_pars$"mi-Eo"* preds_test$invtemp)
##Plot predictions versus generated data
ggplot(preds_test, aes(x=mi2, y=est))+geom_point(color="blue") +geom_point(data=test_dat, aes(x=mi_2, y=y_0.3))+xlab("Metabolic Index")+ylab("Fish Density")+ylim(0,20000)

#Add observation error w/ Tweedie dispersion
preds_test$pred2 <- rTweedie(preds_test$est, p = test_pars$tweedie_p, phi = test_pars$phi)
ggplot(preds_test, aes(x=mi2, y=pred2))+geom_point(color="blue") +geom_point(data=test_dat, aes(x=mi_0.3, y=y_0.3))+xlab("Metabolic Index")+ylab("Fish Density")

##Calculate difference between predictions and simulated data
preds_test$diff <- preds_test$pred2-test_dat$y_0.3

ggplot(preds_test, aes(x=lon, y=lat))+geom_point(aes(color=log(diff)), size=.5) +scale_color_viridis_c()

###Do for all simulations


##Make lists
Eos <- pars_wide$`mi-Eo`
ps <- pars_wide$tweedie_p
phis <- pars_wide$phi
deltas <- pars_wide$"mi-delta"
s50s <-  pars_wide$"mi-s50"
smaxs <- pars_wide$"mi-smax"
depths <- pars_wide$log_depth_sc
depths2 <- pars_wide$log_depth_sc2


#Remove mesh
data_sims_d <- flatten(data_sims)
data_sims_d <- keep(.x=data_sims_d, .p=is.data.frame)

##Apply to each
preds_mi <- mapply(clean_preds, data_sims_d, fits, Eos, ps, phis, deltas, s50s, smaxs, depths, depths2, SIMPLIFY=FALSE)

###95% confidence interval of predictions and metabolic index threshold effects
##Create mesh from original data
mesh_pred <-  make_mesh(pred_dat, xy_cols = c("lon", "lat"), cutoff = .015)
##Prediction from mean, +sd and -sd
pred_mean <- sdmTMB_simulate(formula=~1+logistic(mi)+depth_scaled+depth_scaled2,
                             data=pred_dat,
                             family=tweedie(link="log"),
                             time="year",
                             tweedie_p=pars_performance[11,2],
                             phi=pars_performance[8,2],
                             range=pars_performance[9,2],
                             sigma_O=pars_performance[10,2],
                             sigma_E=NULL,
                             mesh=mesh_pred,
                             threshold_coefs=c(pars_performance[6,2], pars_performance[4,2],pars_performance[7,2], pars_performance[5,2]),
                             B=c(pars_performance[1,2], pars_performance[2,2], pars_performance[3,2]))

pred_low <- sdmTMB_simulate(formula=~1+logistic(mi)+depth_scaled+depth_scaled2,
                            data=pred_dat,
                            family=tweedie(link="log"),
                            time="year",
                            tweedie_p=pars_performance[11,2],
                            phi=pars_performance[8,2],
                            range=pars_performance[9,2],
                            sigma_O=pars_performance[10,2],
                            sigma_E=NULL,
                            mesh=mesh_pred,
                            threshold_coefs=c((pars_performance[6,2]-pars_performance[6,3]), (pars_performance[4,2]-pars_performance[4,3]),(pars_performance[7,2]-pars_performance[7,3]), (pars_performance[5,2]-pars_performance[5,3])),
                            B=c((pars_performance[1,2]-pars_performance[1,2]), (pars_performance[2,2]-pars_performance[2,3]), (pars_performance[3,2]-pars_performance[3,3])))

pred_high <- sdmTMB_simulate(formula=~1+logistic(mi)+depth_scaled+depth_scaled2,
                             data=pred_dat,
                             family=tweedie(link="log"),
                             time="year",
                             tweedie_p=pars_performance[11,2],
                             phi=pars_performance[8,2],
                             range=pars_performance[9,2],
                             sigma_O=pars_performance[10,2],
                             sigma_E=NULL,
                             mesh=mesh_pred,
                             threshold_coefs=c((pars_performance[6,2]+pars_performance[6,3]), (pars_performance[4,2]+pars_performance[4,3]),(pars_performance[7,2]+pars_performance[7,3]), (pars_performance[5,2]+pars_performance[5,3])),
                             B=c((pars_performance[1,2]+pars_performance[1,2]), (pars_performance[2,2]+pars_performance[2,3]), (pars_performance[3,2]+pars_performance[3,3])))

##Combine
pred_int <- as.data.frame(pred_mean$observed)
colnames(pred_int) <- "mean"
pred_int$high <- pred_high$observed
pred_int$low <- pred_low$observed
pred_int$mi <- pred_dat$po2*exp(pars_performance[5,2]* pred_dat$invtemp)
pred_int$mi_low <- pred_dat$po2*exp(pars_performance[5,2]-pars_performance[5,3]* pred_dat$invtemp)
pred_int$mi_high <- pred_dat$po2*exp(pars_performance[5,2]+pars_performance[5,3]* pred_dat$invtemp)

ggplot(pred_int, aes(y=mean,x=mi))+geom_point()+geom_point(aes(y=low,x=mi_low), color="green")+geom_point(aes(y=high,x=mi_high), color="red")

##Comparison of metabolic index effect
pred_int$mi_effect_mean <- (pars_performance[7,2]) * (1 / (1 + exp(-log(19) * (pred_int$mi - pars_performance[6,2]) / pars_performance[4,2])) - 1)
pred_int$mi_effect_low <- (pars_performance[7,2]-pars_performance[7,3]) * (1 / (1 + exp(-log(19) * (pred_int$mi_low - (pars_performance[6,2]-pars_performance[6,3])) / (pars_performance[4,2]-pars_performance[4,3]))) - 1)
pred_int$mi_effect_high <- (pars_performance[7,2]+pars_performance[7,3]) * (1 / (1 + exp(-log(19) * (pred_int$mi_high - (pars_performance[6,2]+pars_performance[6,3])) / (pars_performance[4,2]+pars_performance[4,3]))) - 1)

ggplot(pred_int, aes(x=mi, y=mi_effect_mean))+geom_line()+geom_line(aes(x=mi_low, y=mi_effect_low), color="green")+geom_line(aes(x=mi_high, y=mi_effect_high),color='red')+xlab("Metabolic Index")+ylab("Effect")

#Plot depth vs metabolic index effect
pred_int$depth_scaled <- dat$log_depth_scaled
pred_high$depth_scaled <- dat$log_depth_scaled
pred_low$depth_scaled <- dat$log_depth_scaled
ggplot(pred_int, aes(x=depth_scaled, y=mi_effect_mean))+geom_point(aes(color=log(mean)))+xlab("Depth")+ylab("Metabolic Index Effect")

#Plot metabolic index against depth
ggplot(data=subset(pred_int, mean > 0), aes(x=depth_scaled, y=mi_effect_mean))+geom_point(aes(color=log(mean)))+xlab("Depth")+ylab("Metabolic Index Effect")
#Against depth & depth2 effect
pred_int$depth_effect <- pred_mean$depth_scaled
pred_int$depth2_effect <- pred_mean$depth_scaled2
pred_int$depth_scaled2 <- dat$log_depth_scaled2
ggplot(data=subset(pred_int, mean > 0), aes(x=depth_scaled, y=mi_effect_mean))+geom_point(aes(color=log(mean)))+xlab("Depth Effect")+ylab("Metabolic Index Effect")

ggplot(pred_int, aes(x=mi, y=mi_effect_mean))+geom_line()+geom_point(aes(x=mi, y=depth_effect, color=depth_scaled))+xlab("Metabolic Index")+ylab("Metabolic Index Effect")
ggplot(pred_int, aes(x=depth_scaled, y=mi_effect_mean))+geom_line()+geom_point(aes(x=mi, y=depth_effect, color=depth_scaled))+xlab("Metabolic Index")+ylab("Metabolic Index Effect")

pred_int$depth_effect_combined <- pred_int$depth2_effect+pred_int$depth_effect

pred_int$mi_scaled <- scale(pred_int$mi)
ggplot(pred_int, aes(x=depth_scaled, y=depth_effect_combined))+geom_line()+geom_line(data=pred_int, aes(x=mi_scaled, y=mi_effect_mean), color="purple")+xlab("Data Scaled")+ylab("Effect Size")

ggplot(pred_int, aes(x=depth_effect_combined, y=mi_effect_mean))+geom_point()+xlab("Combined Depth Effect")+ylab("Metabolic Index Effect")

+geom_point(aes(x=mi, y=depth_effect_combined, color=depth_scaled))+xlab("Metabolic Index")+ylab("Metabolic Index Effect")


##With one example
preds_mi_1 <- as.data.frame(preds_mi[20])

#Plot depth vs metabolic index effect
ggplot(preds_mi_1, aes(x=log_depth_sc, y=depth_effect_combined))+geom_point()

#Plot metabolic index against depth
ggplot(preds_mi_1, aes(x=log_depth_sc, y=mi_effect))+geom_point()

ggplot(preds_mi_1, aes(x=depth_effect_combined, y=mi_effect))+geom_point()





###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##For better-fitting models
make_pars_wide <- function(pars){
  pars <- pars[1:3]
  pars <- pivot_wider(pars, id_cols=id, names_from=term, values_from=estimate)
}

pars_wide <- make_pars_wide(pars1)

predict_sims <- function(x){
  if(!is.character(x)){
    preds <- predict(x, newdata=pred_dat, type="response", return_tmb_object=T)
    return <- preds
  }
} 

clean_preds <- function(dat, fit, Eo, p, phi, deltas, s50s, smaxs, depths,depths2){
  pred_dat2 <- dat[,c('po2', 'invtemp', 'log_depth_sc', 'log_depth_sc2', 'lat', 'lon', 'year')]
  preds_test <- predict(fit, newdata=pred_dat2, type="response")
  #preds_test$pred <- exp(preds_test$est)
  ##Calculate metabolic index
  preds_test$mi2 <- preds_test$po2*exp(Eo* preds_test$invtemp)
  ##Plot predictions versus generated data
  #ggplot(preds_test, aes(x=mi2, y=pred))+geom_point(color="blue") +geom_point(data=test_dat, aes(x=mi1, y=y_0.3))+xlab("Metabolic Index")+ylab("Fish Density")
  
  #Add observation error w/ Tweedie dispersion
  preds_test$pred2 <- rTweedie(preds_test$est, p = p, phi = phi)
  preds_test$mi_effect <- (smaxs) * (1 / (1 + exp(-log(19) * (preds_test$mi2 - s50s) / deltas)) - 1)
  preds_test$depth_effect_combined <- (preds_test$log_depth_sc*depths)+(preds_test$log_depth_sc2*depths2)
  return(preds_test)
}

calculate_mi_quant<- function(dat,quant, type){
  if (type==1) {
    total <- sum(dat$y_0.3)
    test_dat <- arrange(dat, mi_0.3)
    test_dat$csum <- cumsum(test_dat$y_0.3)
    test_dat$csum_prop <- test_dat$csum/total
    crit <- test_dat[which.min(abs(quant-test_dat$csum_prop)),]
    crit <- crit$mi_0.3
  }
  if (type==2) {
    total <- sum(dat$pred2)
    test_dat <- arrange(dat, mi2)
    test_dat$csum <- cumsum(test_dat$pred2)
    test_dat$csum_prop <- test_dat$csum/total
    crit <- test_dat[which.min(abs(quant-test_dat$csum_prop)),]$mi2
  }
  if (type==3) {
    total <- sum(dat$pred)
    test_dat <- arrange(dat, mi2)
    test_dat$csum <- cumsum(test_dat$pred)
    test_dat$csum_prop <- test_dat$csum/total
    crit <- test_dat[which.min(abs(quant-test_dat$csum_prop)),]$mi2
  }
  return(crit)
}

#####~~~~~~~~~~~~Model Predictions~~~~~~~~~
###Make dataset to make predictions on (original true sablefish data)
kelvin = 273.15
boltz = 0.000086173324
tref <- 12
pred_dat <- dat[,c('po2', 'temp', 'depth', 'X', 'Y', 'year')]
pred_dat$lon <- pred_dat$X/1000
pred_dat$lat <- pred_dat$Y/1000
pred_dat$log_depth_sc <- scale(log(pred_dat$depth))
pred_dat$log_depth_sc2 <- pred_dat$log_depth_sc ^ 2
pred_dat$depth <- NULL
pred_dat$invtemp <- (1 / boltz)  * ( 1 / (pred_dat$temp + 273.15) - 1 / (tref + 273.15))
pred_dat$temp <- NULL
pred_dat$latitude <- NULL
pred_dat$longitude <- NULL

###Apply to all model fits
preds<- lapply(fits, predict_sims)

###Example of just one
##Prediction
test2 <- preds[[20]]
test <- test2[[1]]

##Predict from simulated data and model fit
test_fit <- fits[[20]]
##Simulated data
test_dat <- data_sims[[20]]
test_dat <- as.data.frame(test_dat[1])
##Estimated parameters
test_pars <- pars_wide[20,]
##Calculate metabolic index
test_dat$mi_2 <- test_dat$po2*exp(test_pars$"mi-Eo"* test_dat$invtemp)

##Prediction from simulated data, not real data
pred_dat2 <- test_dat[,c('po2', 'invtemp', 'log_depth_sc', 'log_depth_sc2', 'lat', 'lon', 'year')]
preds_test <- predict(test_fit, newdata=pred_dat2, type="response", return_tmb_object=T)
preds_test <- preds_test[[1]]
preds_test$pred <- exp(preds_test$est)
##Calculate metabolic index
preds_test$mi2 <- preds_test$po2*exp(test_pars$"mi-Eo"* preds_test$invtemp)
##Plot predictions versus generated data
ggplot(preds_test, aes(x=mi2, y=est))+geom_point(color="blue") +geom_point(data=test_dat, aes(x=mi_2, y=y_0.3))+xlab("Metabolic Index")+ylab("Fish Density")+ylim(0,20000)

#Add observation error w/ Tweedie dispersion
preds_test$pred2 <- rTweedie(preds_test$est, p = test_pars$tweedie_p, phi = test_pars$phi)
ggplot(preds_test, aes(x=mi2, y=pred2))+geom_point(color="blue") +geom_point(data=test_dat, aes(x=mi_0.3, y=y_0.3))+xlab("Metabolic Index")+ylab("Fish Density")

##Calculate difference between predictions and simulated data
preds_test$diff <- preds_test$pred2-test_dat$y_0.3

ggplot(preds_test, aes(x=lon, y=lat))+geom_point(aes(color=log(diff)), size=.5) +scale_color_viridis_c()

###Do for all simulations


##Make lists
Eos <- pars_wide$`mi-Eo`
ps <- pars_wide$tweedie_p
phis <- pars_wide$phi
deltas <- pars_wide$"mi-delta"
s50s <-  pars_wide$"mi-s50"
smaxs <- pars_wide$"mi-smax"
depths <- pars_wide$log_depth_sc
depths2 <- pars_wide$log_depth_sc2


#Remove mesh
data_sims_d <- flatten(data_sims)
data_sims_d <- keep(.x=data_sims_d, .p=is.data.frame)

##Apply to each
preds_mi <- mapply(clean_preds, data_sims_d, fits, Eos, ps, phis, deltas, s50s, smaxs, depths, depths2, SIMPLIFY=FALSE)

###95% confidence interval of predictions and metabolic index threshold effects
##Create mesh from original data
mesh_pred <-  make_mesh(pred_dat, xy_cols = c("lon", "lat"), cutoff = .015)
##Prediction from mean, +sd and -sd
pred_mean <- sdmTMB_simulate(formula=~1+logistic(mi)+depth_scaled+depth_scaled2,
                             data=pred_dat,
                             family=tweedie(link="log"),
                             time="year",
                             tweedie_p=pars_performance[11,2],
                             phi=pars_performance[8,2],
                             range=pars_performance[9,2],
                             sigma_O=pars_performance[10,2],
                             sigma_E=NULL,
                             mesh=mesh_pred,
                             threshold_coefs=c(pars_performance[6,2], pars_performance[4,2],pars_performance[7,2], pars_performance[5,2]),
                             B=c(pars_performance[1,2], pars_performance[2,2], pars_performance[3,2]))

pred_low <- sdmTMB_simulate(formula=~1+logistic(mi)+depth_scaled+depth_scaled2,
                            data=pred_dat,
                            family=tweedie(link="log"),
                            time="year",
                            tweedie_p=pars_performance[11,2],
                            phi=pars_performance[8,2],
                            range=pars_performance[9,2],
                            sigma_O=pars_performance[10,2],
                            sigma_E=NULL,
                            mesh=mesh_pred,
                            threshold_coefs=c((pars_performance[6,2]-pars_performance[6,3]), (pars_performance[4,2]-pars_performance[4,3]),(pars_performance[7,2]-pars_performance[7,3]), (pars_performance[5,2]-pars_performance[5,3])),
                            B=c((pars_performance[1,2]-pars_performance[1,2]), (pars_performance[2,2]-pars_performance[2,3]), (pars_performance[3,2]-pars_performance[3,3])))

pred_high <- sdmTMB_simulate(formula=~1+logistic(mi)+depth_scaled+depth_scaled2,
                             data=pred_dat,
                             family=tweedie(link="log"),
                             time="year",
                             tweedie_p=pars_performance[11,2],
                             phi=pars_performance[8,2],
                             range=pars_performance[9,2],
                             sigma_O=pars_performance[10,2],
                             sigma_E=NULL,
                             mesh=mesh_pred,
                             threshold_coefs=c((pars_performance[6,2]+pars_performance[6,3]), (pars_performance[4,2]+pars_performance[4,3]),(pars_performance[7,2]+pars_performance[7,3]), (pars_performance[5,2]+pars_performance[5,3])),
                             B=c((pars_performance[1,2]+pars_performance[1,2]), (pars_performance[2,2]+pars_performance[2,3]), (pars_performance[3,2]+pars_performance[3,3])))

##Combine
pred_int <- as.data.frame(pred_mean$observed)
colnames(pred_int) <- "mean"
pred_int$high <- pred_high$observed
pred_int$low <- pred_low$observed
pred_int$mi <- pred_dat$po2*exp(pars_performance[5,2]* pred_dat$invtemp)
pred_int$mi_low <- pred_dat$po2*exp(pars_performance[5,2]-pars_performance[5,3]* pred_dat$invtemp)
pred_int$mi_high <- pred_dat$po2*exp(pars_performance[5,2]+pars_performance[5,3]* pred_dat$invtemp)

ggplot(pred_int, aes(y=mean,x=mi))+geom_point()+geom_point(aes(y=low,x=mi_low), color="green")+geom_point(aes(y=high,x=mi_high), color="red")

##Comparison of metabolic index effect
pred_int$mi_effect_mean <- (pars_performance[7,2]) * (1 / (1 + exp(-log(19) * (pred_int$mi - pars_performance[6,2]) / pars_performance[4,2])) - 1)
pred_int$mi_effect_low <- (pars_performance[7,2]-pars_performance[7,3]) * (1 / (1 + exp(-log(19) * (pred_int$mi_low - (pars_performance[6,2]-pars_performance[6,3])) / (pars_performance[4,2]-pars_performance[4,3]))) - 1)
pred_int$mi_effect_high <- (pars_performance[7,2]+pars_performance[7,3]) * (1 / (1 + exp(-log(19) * (pred_int$mi_high - (pars_performance[6,2]+pars_performance[6,3])) / (pars_performance[4,2]+pars_performance[4,3]))) - 1)

ggplot(pred_int, aes(x=mi, y=mi_effect_mean))+geom_line()+geom_line(aes(x=mi_low, y=mi_effect_low), color="green")+geom_line(aes(x=mi_high, y=mi_effect_high),color='red')+xlab("Metabolic Index")+ylab("Effect")

#Plot depth vs metabolic index effect
pred_int$depth_scaled <- dat$log_depth_scaled
pred_high$depth_scaled <- dat$log_depth_scaled
pred_low$depth_scaled <- dat$log_depth_scaled
ggplot(pred_int, aes(x=depth_scaled, y=mi_effect_mean))+geom_point(aes(color=log(mean)))+xlab("Depth")+ylab("Metabolic Index Effect")

#Plot metabolic index against depth
ggplot(data=subset(pred_int, mean > 0), aes(x=depth_scaled, y=mi_effect_mean))+geom_point(aes(color=log(mean)))+xlab("Depth")+ylab("Metabolic Index Effect")
#Against depth & depth2 effect
pred_int$depth_effect <- pred_mean$depth_scaled
pred_int$depth2_effect <- pred_mean$depth_scaled2
pred_int$depth_scaled2 <- dat$log_depth_scaled2
ggplot(data=subset(pred_int, mean > 0), aes(x=depth_scaled, y=mi_effect_mean))+geom_point(aes(color=log(mean)))+xlab("Depth Effect")+ylab("Metabolic Index Effect")

ggplot(pred_int, aes(x=mi, y=mi_effect_mean))+geom_line()+geom_point(aes(x=mi, y=depth_effect, color=depth_scaled))+xlab("Metabolic Index")+ylab("Metabolic Index Effect")
ggplot(pred_int, aes(x=depth_scaled, y=mi_effect_mean))+geom_line()+geom_point(aes(x=mi, y=depth_effect, color=depth_scaled))+xlab("Metabolic Index")+ylab("Metabolic Index Effect")

pred_int$depth_effect_combined <- pred_int$depth2_effect+pred_int$depth_effect

pred_int$mi_scaled <- scale(pred_int$mi)
ggplot(pred_int, aes(x=depth_scaled, y=depth_effect_combined))+geom_line()+geom_line(data=pred_int, aes(x=mi_scaled, y=mi_effect_mean), color="purple")+xlab("Data Scaled")+ylab("Effect Size")

ggplot(pred_int, aes(x=depth_effect_combined, y=mi_effect_mean))+geom_point()+xlab("Combined Depth Effect")+ylab("Metabolic Index Effect")

+geom_point(aes(x=mi, y=depth_effect_combined, color=depth_scaled))+xlab("Metabolic Index")+ylab("Metabolic Index Effect")


##With one example
preds_mi_1 <- as.data.frame(preds_mi[20])

#Plot depth vs metabolic index effect
ggplot(preds_mi_1, aes(x=log_depth_sc, y=depth_effect_combined))+geom_point()

#Plot metabolic index against depth
ggplot(preds_mi_1, aes(x=log_depth_sc, y=mi_effect))+geom_point()

ggplot(preds_mi_1, aes(x=depth_effect_combined, y=mi_effect))+geom_point()


