###Install packages
remotes::install_github("pbs-assess/sdmTMB", dependencies = TRUE, ref="mi", force=TRUE)
library(sdmTMB)
library(ggplot2)
library(visreg)
library(ggeffects)
library(dplyr)
library(tidyr)
library(mgcv)

###Data
here("thresholds_mi_distribution")
setwd("~/Dropbox/Mac/Documents/GitHub/estimating_mi_from_distribution")
dat <- readRDS("data_sablefish2.rds")

##Re-name MI just to be safe with the sdmTMB modeling
dat$metabolic_index <- dat$mi
dat$mi <- NULL

##Calculate inverse temp
kelvin = 273.15 #To convert to Kelvin
boltz = 0.000086173324 #Boltzman's constant
tref <- 12 #Reference temperature in celsius
#Calculate inverse temp
dat$invtemp <- (1 / boltz)  * ( 1 / (dat$temp + 273.15) - 1 / (tref + 273.15))

###Set ggplot theme
theme_set(theme_bw(base_size = 15))
theme_update(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

###Look at data
ggplot(dat, aes(x=temp, y=po2))+geom_point()

###Make mesh
mesh <-
  make_mesh(data = dat,
            xy_cols = c("X", "Y"),
            n_knots = 250) # choose # knots

###Alternative models

##Set up delta-AIC table
models <- c("breakpt(o2)", "internal Eo prior", 'internal Eo no prior', "logistic(po2)", "logistic(mi)", "breakpt(temp-corrected o2 from prior", "breakpt(temp-corrected o2 from no prior")
AIC <- as.data.frame(matrix(NA, ncol = 1, nrow = 7, dimnames = list(models)))


##Breakpoint pO2
m1 <- sdmTMB(cpue_kg_km2 ~ -1+as.factor(year)+breakpt(po2_sc)+log_depth_scaled+log_depth_scaled2, 
             data = dat, 
             spatial = "on",
             spatiotemporal="off",
             time=NULL,
             reml=F,
             mesh=mesh,
             anisotropy=T,
             family =tweedie(link="log"),
             control = sdmTMBcontrol(
               nlminb_loops = 2))

summary(m1)
AIC[1,] <- AIC(m1)

##Internal estimation of Eo and logistic on temp-corrected o2 w/ prior on Eo
#Set starting parameters: 
#Note: have to start it at least above 0 (e.g. 0.05) or results in NaN gradient evaluation, but starting
#at 2 vs 0.05, for instance, doesn't change parameter estimation or maximum final gradient
start <- matrix(0, ncol = 1, nrow = 4)
start[1, 1] <- 1 #s50
start[2, 1] <- 1 #delta
start[3, 1] <- 10 #smax (ie beta_3)
start[4, 1] <- 0.448 #Eo--this is the value used in the paper
m2 <- sdmTMB(cpue_kg_km2 ~ 1+year+logistic(mi)+log_depth_scaled+log_depth_scaled2, 
             data = dat, 
             spatial = "on",
             mesh=mesh,
             family =tweedie(link="log"),
             priors=sdmTMBpriors(threshold = normal(c(NA, NA, 10, 0.448), c(NA, NA, 1, 0.15))),
             control = sdmTMBcontrol(
               start = list(b_threshold = start),
               lower = list(b_threshold = c(0, 0, 0, 0)), upper = list(b_threshold = c(Inf, Inf,50, Inf)),
               newton_loops = 1))

summary(m2)
AIC(m2)
AIC[2,] <- AIC(m2)

##Internal estimation of Eo and logistic on temp-corrected o2 without prior
start <- matrix(0, ncol = 1, nrow = 4)
start[1, 1] <- 1 #s50
start[2, 1] <- 1 #delta
start[3, 1] <- 10 #smax (ie beta_3)
start[4, 1] <- 0.448 #Eo
m3 <- sdmTMB(cpue_kg_km2 ~ -1+year+logistic(mi)+log_depth_scaled+log_depth_scaled2, 
             data = dat, 
             spatial = "on",
             mesh=mesh,
             family =tweedie(link="log"),
             anisotropy=T,
             reml=F,
             time=NULL,
             priors=sdmTMBpriors(threshold = normal(c(NA, NA, NA, NA), c(NA, NA, NA, NA))),
             control = sdmTMBcontrol(
               start = list(b_threshold = start),
               lower = list(b_threshold = c(0, 0, 0, 0)), upper = list(b_threshold = c(Inf, Inf,50, Inf)),
               newton_loops = 1))

summary(m3)
AIC(m3)
AIC[3,] <- AIC(m3)

##Logistic/sigmoid on po2
start <- matrix(0, ncol = 1, nrow = 3)
start[1, 1] <- -5 #s50
start[2, 1] <- -0.3 #s95
start[3, 1] <- 1 #smax
m4 <- sdmTMB(cpue_kg_km2 ~ -1+year+logistic(po2_sc)+log_depth_scaled+log_depth_scaled2, 
             data = dat, 
             spatial = "on",
             mesh=mesh,
             anisotropy=T,
             reml=F,
             time=NULL,
             family =tweedie(link="log"),
             control = sdmTMBcontrol(
              start = list(b_threshold=start),
               newton_loops = 2,
             eval.max = 1000000L,
             iter.max = 1000000L,
             nlminb_loops = 100L,
             lower = list(b_threshold = c(-6, -1, 0)),
             upper= list(b_threshold=c(Inf, Inf, Inf))
             ))
     
summary(m4)
AIC(m4)
AIC[4,] <- AIC(m4)

##Logistic/sigmoid on metabolic index (from the values used in the paper)
start <- matrix(0, ncol = 1, nrow = 3)
start[1, 1] <- 0.1 #s50
start[2, 1] <- 0.3#s95
start[3, 1] <- 0.1 #smax
m5 <- sdmTMB(cpue_kg_km2 ~ -1+year+logistic(metabolic_index)+log_depth_scaled+log_depth_scaled2, 
             data = dat, 
             spatial = "on",
             mesh=mesh,
             anisotropy=T,
             reml=F,
             time=NULL,
             family =tweedie(link="log"),
             control = sdmTMBcontrol(
               newton_loops = 1,
             eval.max = 1000000L,
             iter.max = 1000000L,
             nlminb_loops = 100L,
             lower = list(b_threshold = c(0.1, 0.2, 0)),
             upper= list(b_threshold=c(Inf, Inf, Inf))
))

summary(m5)
AIC(m5)
AIC[5,] <- AIC(m5)

##Calculate temperature-corrected oxygen w/ Eo value estimated from the model
#Extract Eo
Eo_m2 <- as.numeric(m2[["sd_report"]][["par.fixed"]][16])
Eo_m3 <- as.numeric(m3[["sd_report"]][["par.fixed"]][16])

#Calculate
dat$metabolic_index_est_m2 = dat$po2*exp(Eo_m2 * dat$invtemp)
dat$metabolic_index_est_m3 = dat$po2*exp(Eo_m3 * dat$invtemp)
dat$metabolic_index_est_m2 <- scale(dat$metabolic_index_est_m2)
dat$metabolic_index_est_m3 <- scale(dat$metabolic_index_est_m3)

#Breakpoint models with calculated MI
m6 <- sdmTMB(cpue_kg_km2 ~ -1+year+breakpt(metabolic_index_est_m2)+log_depth_scaled+log_depth_scaled2, 
             data = dat, 
             spatial = "on",
             mesh=mesh,
             anisotropy=T,
             reml=F,
             time=NULL,
             family =tweedie(link="log"),
             control = sdmTMBcontrol(
               newton_loops = 1))

summary(m6)
AIC(m6)
AIC[6,] <- AIC(m6)

m7 <- sdmTMB(cpue_kg_km2 ~ -1+year+breakpt(metabolic_index_est_m3)+log_depth_scaled+log_depth_scaled2, 
             data = dat, 
             spatial = "on",
             mesh=mesh,
             family =tweedie(link="log"),
             control = sdmTMBcontrol(
               newton_loops = 1))

summary(m7)
AIC(m7)
AIC[7,] <- AIC(m7)

##Calculate delta-AIC
AIC$dAIC <- abs(min(AIC$V1)-(AIC$V1))

###Plot marginal effects
##Best-fitting internal Eo model 
#Extract parameters
par_estimates <- as.data.frame(tidy(m3, conf.int = TRUE, effects="fixed"))
par_estimates_rand <- as.data.frame(tidy(m3, conf.int = TRUE, effects="ran_pars"))
par_estimates <- bind_rows(par_estimates, par_estimates_rand)
pars <- par_estimates[1:2]
pars <- pivot_wider(pars, names_from=term, values_from=estimate)
#Calculate temp-corrected po2 from Eo value
dat$mi_pred <- dat$po2*exp(pars$"mi-Eo"* dat$invtemp)
#Calculate effect of MI
dat$mi_effect <- pars$"mi-smax" * (1 / (1 + exp(-log(19) * (dat$mi_pred - pars$`mi-s50`) / pars$"mi-delta")) - 1)
#Calculate combined depth effect
dat$depth_effect_combined <- (dat$log_depth_scaled*pars$log_depth_scaled)+(dat$log_depth_scaled2*pars$log_depth_scaled2)

##For breakpt(o2) model
par_estimates2 <- as.data.frame(tidy(m1, conf.int = TRUE, effects="fixed"))
par_estimates_rand2 <- as.data.frame(tidy(m1, conf.int = TRUE, effects="ran_pars"))
par_estimates2 <- bind_rows(par_estimates2, par_estimates_rand2)
pars2 <- par_estimates2[1:2]
pars2 <- pivot_wider(pars2, names_from=term, values_from=estimate)
#Calculate depth effects
dat$depth_effect_combined2 <- (dat$log_depth_scaled*pars2$log_depth_scaled)+(dat$log_depth_scaled2*pars2$log_depth_scaled2)
#Calculate po2 effect
dat$po2_effect <- ifelse(dat$po2_sc>pars2$`po2_sc-breakpt`,pars2$`po2_sc-breakpt`*pars2$`po2_sc-slope`,dat$`po2_sc`*pars2$`po2_sc-slope`)

#Plot depth effects compared between model
ggplot(dat, aes(y=depth_effect_combined, x=depth))+geom_point()+geom_point(dat, mapping=aes(y=depth_effect_combined2, x=depth), color="red")

#Plot MI effect (for internal Eo model)
ggplot(dat, aes(x=metabolic_index, y=mi_effect))+geom_point()+geom_point(mapping=aes(x=metabolic_index, y=depth_effect_combined))
#Color with depth
ggplot(subset(dat, metabolic_index<3), aes(x=metabolic_index, y=mi_effect))+geom_point()
#subset to low MI
ggplot(subset(dat, metabolic_index<3), aes(x=metabolic_index, y=mi_effect))+geom_point()

#Plot breakpt(po2) effect
ggplot(dat, aes(x=po2, y=po2_effect))+geom_point()

##Other figures for fun--these can be ignored
#Plot depth vs metabolic index effect
ggplot(dat, aes(x=log_depth_scaled, y=depth_effect_combined))+geom_point()
#Plot metabolic index against depth
ggplot(dat, aes(x=log_depth_scaled, y=mi_effect))+geom_point()
ggplot(dat, aes(x=depth_effect_combined, y=mi_effect))+geom_point()
#Subset MI 0-5/7
dat_sub <- subset(dat, mi_pred<4)
ggplot(dat_sub, aes(x=log_depth_scaled, y=mi_effect))+geom_point()
ggplot(dat_sub, aes(x=mi_pred, y=mi_effect))+geom_point()
ggplot(dat_sub, aes(x=depth_effect_combined, y=mi_effect))+geom_point()
ggplot(dat_sub, aes(x=mi_pred, y=depth_effect_combined))+geom_point()
#Plot metabolic index effect against depth
ggplot(dat, aes(x=log_depth_scaled, y=mi_effect))+geom_point()
ggplot(dat, aes(x=depth_effect_combined, y=mi_effect))+geom_point()

###DHARMa residuals
#Predictions (because DHARMa requires them)
preds2 <- predict(m3, newdata=dat, type="response", return_tmb_object=T)
preds <-  as.data.frame(preds2[1])

#plot against real data
ggplot(preds, aes(x=data.mi_pred, y=data.est))+geom_point(color="blue") +geom_point(data=dat, aes(x=mi_pred, y=cpue_kg_km2))+xlab("Metabolic Index")+ylab("Fish Density")+ylim(0,20000)

#dharma residuals
s_res <- simulate(m3, nsim = 500)
pred_fixed <- m3$family$linkinv(predict(m3)$est_non_rf)
dharma <- DHARMa::createDHARMa(
  simulatedResponse = s_res,
  observedResponse = dat$cpue_kg_km2,
  fittedPredictedResponse = pred_fixed
)

plot(dharma)
DHARMa::testResiduals(dharma)
resids <- dharma$scaledResiduals
dat$resids <- resids

#Plot residuals
ggplot(dat, (aes(x=log_depth_scaled, y=resids)))+geom_point()
ggplot(dat, (aes(x=log_depth_scaled2, y=resids)))+geom_point()

#Plot spatial random variation
ggplot(preds, aes(data.X, data.Y, col = data.omega_s)) +
  scale_colour_gradient2() +
  geom_point(size=0.5)

#Other stuff with predictions
#Add in observation error
preds$pred2 <- rTweedie(preds$data.est, p = pars$tweedie_p, phi = pars$phi)
#Plot spatial variation effect
ggplot(preds, aes(x=data.omega_s, y=log_depth_scaled2))+geom_point()

















########More versions fit in early March--ignore this stuff##############
###Fit models
#Breakpoint function on pre-calculated MI
m1 <- sdmTMB(cpue_kg_km2 ~ 1+year+breakpt(metabolic_index_sc)+log_depth_scaled+log_depth_scaled2, 
            data = dat, 
            spatial = "on",
            mesh=mesh,
            family =tweedie(link="log"),
            control=sdmTMBcontrol(newton_loops = 1))
#If use raw metabolic index: non-positive definite hessian, Note: 0 and 0 estimates for breakpt slope and threshold
summary(m1)
AIC(m1)

#Breakpoint function on curved Eo
m2 <- sdmTMB(cpue_kg_km2 ~ 1+year+breakpt(mi_curve_sc)+log_depth_scaled+log_depth_scaled2, 
             data = dat, 
             spatial = "on",
             mesh=mesh,
             family =tweedie(link="log"),
             control = sdmTMBcontrol(
               newton_loops = 1))
#Note: if use raw mi_curve: non-positive definite hessian, 0 and 0 estimates for breakpt slope and threshold
summary(m2)
AIC(m2)

##Estimate Eo and MI
#Set starting parameters: 
#Note: have to start it at least above 0 (e.g. 0.05) or results in NaN gradient evaluation, but starting
#at 2 vs 0.05, for instance, doesn't change parameter estimation or maximum final gradient
#
start <- matrix(0, ncol = 1, nrow = 4)
start[1, 1] <- 0.05 #s50
start[2, 1] <- 0.05 #delta
start[3, 1] <- 1 #smax (ie beta_3)
start[4, 1] <- 0.3 #Eo
m3 <- sdmTMB(cpue_kg_km2 ~ 1+year+logistic(mi)+log_depth_scaled+log_depth_scaled2, 
             data = dat, 
             spatial = "on",
             mesh=mesh,
             family =tweedie(link="log"),
             priors=sdmTMBpriors(threshold = normal(c(NA, NA, NA, NA), c(NA, NA, NA, NA))),
             control = sdmTMBcontrol(
               start = list(b_threshold = start),
               lower = list(b_threshold = c(0, 0, 0, 0)), upper = list(b_threshold = c(Inf, Inf,50, Inf)),
               newton_loops = 1))
             
summary(m3)
AIC(m3)

#Breakpt on o2, no metabolic index
m5 <- sdmTMB(cpue_kg_km2 ~ 1+year+breakpt(po2_sc)+log_depth_scaled+log_depth_scaled2, 
             data = dat, 
             spatial = "on",
             mesh=mesh,
             family =tweedie(link="log"),
             control = sdmTMBcontrol(
               newton_loops = 1,
            priors))

summary(m5)
AIC(m5)

start <- matrix(0, ncol = 1, nrow = 4)
start[1, 1] <- 0.05 #s50
start[2, 1] <- 0.05 #delta
start[3, 1] <-5 #smax (ie beta_3)
start[4, 1] <- 0.05 #Eo
m6 <- sdmTMB(cpue_kg_km2 ~ 1+year+logistic(mi)+s(log_depth_scaled), 
             data = dat, 
             spatial = "on",
             mesh=mesh,
             family =tweedie(link="log"),
             priors=sdmTMBpriors(threshold = normal(c(1, 1, 20, 0.3), c(0.25, 0.25, 1, 0.1))),
             control = sdmTMBcontrol(
               start = list(b_threshold = start),
               lower = list(b_threshold = c(0, 0, 0, 0)), upper = list(b_threshold = c(Inf, Inf, 50, Inf)),
               newton_loops = 1))

summary(m6)
AIC(m6)
visreg(m6, "log_depth_scaled")

m7 <- sdmTMB(cpue_kg_km2 ~ 1+year+s(mi_curve_sc)+log_depth_scaled+log_depth_scaled2, 
             data = dat, 
             spatial = "on",
             mesh=mesh,
             family =tweedie(link="log"),
             control = sdmTMBcontrol(
               newton_loops = 1))
summary(m7)
AIC(m7)

visreg(m6, "log_depth_scaled")
visreg(m7, "mi_curve_sc")

m8 <- sdmTMB(cpue_kg_km2 ~ 1+year+s(metabolic_index_sc)+log_depth_scaled+log_depth_scaled2, 
             data = dat, 
             spatial = "on",
             mesh=mesh,
             family =tweedie(link="log"),
             control = sdmTMBcontrol(
               newton_loops = 1))
summary(m8)
AIC(m8)

visreg(m8, "metabolic_index_sc")


