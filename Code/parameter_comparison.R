###~~~~~~~~~~~~~~~~~~Evaluate parameter estimation~~~~~~~~~~~~~~~~~~~~~~~~~
library(ggplot2)
install.packages("dplyr")
library(dplyr)
library(tidyr)
library(ggrepel)
library(broom)
library(purrr)
install.packages("directlabels")
library(directlabels)
install.packages("remotes")
remotes::install_github("teunbrand/ggh4x")
library(ggh4x)
##Set ggplot themes
theme_set(theme_bw(base_size = 15))
theme_update(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#Load models
load("~/Dropbox/Mac/Documents/GitHub/thresholds_mi_distribution/thresholds_mi_distribution/current/output/model_fits.Rdata")

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

###Comparing how starting values impacted parameter estimation for model with no prior
##Function to extract from all fits
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
##Apply to each sim in each model
#Model starting threshold values at correct parameter values
pars1 <- lapply(fits, extract_pars)
#Model starting threshold values lower (at 10% of true value)
pars3 <- lapply(fits3, extract_pars)

##Function to clean up pars for plotting
clean_pars <- function(pars, fits){
  names(pars) <- c(1:length(fits))
  #Remove models with errors
  pars <- keep(pars, function(x) !is.logical(x))
  #Combine into single dataframe, with column of simulation number
  pars <- bind_rows(pars,.id="id")
  return(pars)
}

##Apply to each
#pars1 <- clean_pars(pars1, fits=fits)
#pars2 <- clean_pars(pars3, fits=fits2)

##Add column to label model
#pars1$model <- "Start at true value"
#pars2$model <- "Start at 10% true value"

##Merge together for plotting
#pars <- rbind(pars1,pars2)
#Boxplot of estimates compared to true value
#ggplot(pars, aes(y=estimate, x=term))+geom_boxplot(aes(group=model, fill=model))+facet_wrap("term", scales="free")+theme(legend.position = "none")

###NOTE: Because no difference, going to use fits for rest of comparison
#Make list of fits

fits_all <- list(fits, fits2, fits3)

##Create list of names of models
model_names <- c("Data limited (no prior)", "Data Rich (Correct prior)", "Unusual Case (Incorrect prior)")

##Set true pars vectors of values and names
true_pars2 <- data.frame(term=c("(Intercept)", "log_depth_sc", "log_depth_sc2", "mi-delta", "mi-s50", "mi-smax", "range", "sigma_O", "phi", "tweedie_p", "mi-Eo", "as.factor(year)2011","as.factor(year)2012", "as.factor(year)2013", "as.factor(year)2014", "as.factor(year)2015"), estimate=c(beta0, beta1, beta2, delta, x50, beta3, range, sigma_O, phi, p, 0.3, b_years))
pars_names <- c("(Intercept)","as.factor(year)2011","as.factor(year)2012", "as.factor(year)2013", "as.factor(year)2014", "as.factor(year)2015", "log_depth_sc", "log_depth_sc2", "mi-s50", "mi-delta", "mi-smax", "mi-Eo", "range", "phi", "sigma_O", "tweedie_p") 

##Extract pars for the other two models
##Apply to each sim in each model
pars1 <- lapply(fits, extract_pars)
pars2 <- lapply(fits2, extract_pars)
pars3 <- lapply(fits3, extract_pars)

##Apply to each
pars1 <- clean_pars(pars1, fits=fits)
pars2 <- clean_pars(pars2, fits=fits2)
pars3 <- clean_pars(pars3, fits=fits3)

##Add column to label model
pars1$model <- model_names[1]
pars2$model <- model_names[2]
pars3$model <- model_names[3]

##Merge together for plotting
pars <- rbind(pars1,pars2,  pars3)

##Save
saveRDS(pars, "parameter_estimates.rds")
#Boxplot of estimates compared to true value
ggplot(pars, aes(y=estimate, x=term))+geom_boxplot(aes(group=model, fill=model))+facet_wrap("term", scales="free")+geom_hline(data = true_pars2, aes(yintercept = estimate),linetype="dashed", size=1.2)+theme(legend.position="left", strip.text = element_blank())
#Just Eo
pars$model <- str_wrap(pars$model, width=10)
ggplot(subset(pars, term=='mi-Eo'), aes(y=estimate, x=model))+geom_boxplot(aes(group=model, fill=model))+geom_hline(data = subset(true_pars2, term=="mi-Eo"), aes(yintercept = estimate),linetype="dashed", size=1.2)+theme(legend.position="none")

###Parameter performance measures
##Accuracy
#Average value (accuracy)
avg <- aggregate(estimate ~ term+model, pars, FUN=mean)
#SD of average estimate (precision)
sd_avg <- aggregate(estimate ~ term+model, pars, FUN=stats::sd)
#Range of average estimates (precision)
range_pars <- aggregate(estimate ~ term+model, pars, FUN=fivenum)
#Range of standard error (precision)
sd_range <- aggregate(std.error~ term+model, pars, FUN=fivenum)
#Average standard error (precision)--also to compare to SD of prior
sd_avg <- aggregate(std.error ~ term+model, pars, FUN=mean)
#Standard deviation of standard error 
sd_sd <- aggregate(std.error ~ term+model, pars, FUN=stats::sd)

#Accuracy
#Mean error (model estimate -- truth)
pars$error <- case_when(pars$term=="mi-s50"~pars$estimate-x50,
                        pars$term=="mi-delta"~pars$estimate-delta,
                        pars$term=="mi-smax"~pars$estimate-beta3,
                        pars$term=="mi-Eo"~pars$estimate-Eo,
                        pars$term=="(Intercept)"~pars$estimate-beta0,
                        pars$term=="log_depth_sc"~pars$estimate-beta1,
                        pars$term=="log_depth_sc2"~pars$estimate-beta2,
                        pars$term=="range"~pars$estimate-range,
                        pars$term=="sigma_O"~pars$estimate-sigma_O,
                        pars$term=="phi"~pars$estimate-phi,
                        pars$term=="tweedie_p"~pars$estimate-p,
                        pars$term=="as.factor(year)2011"~pars$estimate-b_years[1],
                        pars$term=="as.factor(year)2012"~pars$estimate-b_years[2],
                        pars$term=="as.factor(year)2013"~pars$estimate-b_years[3],
                        pars$term=="as.factor(year)2014"~pars$estimate-b_years[4],
                        pars$term=="as.factor(year)2015"~pars$estimate-b_years[5])

#Mean squared error
pars$error2 <- pars$error^2
error2 <- aggregate(error2 ~ term+model, pars, FUN=mean)

##Precision
#95% range (upper 95% CI--lower 95% CI)
pars$con.range <- pars$conf.high-pars$conf.low

#Average 95% range
avg_conf.range <- aggregate(con.range ~ term+model, pars, FUN=mean)

##Coverage
#How often does true parameter value fall within the range of the 95% CI
#1 if true value is within confidence interval, 0 if not
pars$within_95 <- case_when(pars$term=="mi-s50"~ifelse(pars$conf.low <x50 & pars$conf.high >x50, 1, 0),
                            pars$term=="mi-delta"~ifelse(pars$conf.low <delta & pars$conf.high >delta, 1, 0),
                            pars$term=="mi-smax"~ifelse(pars$conf.low <beta3 & pars$conf.high >beta3, 1, 0),
                            pars$term=="mi-Eo"~ifelse(pars$conf.low <Eo & pars$conf.high>Eo, 1, 0),
                            pars$term=="(Intercept)"~ifelse(pars$conf.low <beta0 & pars$conf.high>beta0, 1, 0),
                            pars$term=="log_depth_sc"~ifelse(pars$conf.low <beta1 & pars$conf.high>beta1, 1, 0),
                            pars$term=="log_depth_sc2"~ifelse(pars$conf.low <beta2 & pars$conf.high>beta2, 1, 0),
                            pars$term=="range"~ifelse(pars$conf.low <range & pars$conf.high>range, 1, 0),
                            pars$term=="sigma_O"~ifelse(pars$conf.low <sigma_O & pars$conf.high>sigma_O, 1, 0),
                            pars$term=="phi"~ifelse(pars$conf.low <phi & pars$conf.high>phi, 1, 0),
                            pars$term=="tweedie_p"~ifelse(pars$conf.low <p & pars$conf.high>p, 1, 0),
                            pars$term=="as.factor(year)2011"~ifelse(pars$conf.low <b_years[1] & pars$conf.high>b_years[1], 1, 0),
                            pars$term=="as.factor(year)2012"~ifelse(pars$conf.low <b_years[2] & pars$conf.high>b_years[2], 1, 0),
                            pars$term=="as.factor(year)2013"~ifelse(pars$conf.low <b_years[3] & pars$conf.high>b_years[3], 1, 0),
                            pars$term=="as.factor(year)2014"~ifelse(pars$conf.low <b_years[4] & pars$conf.high>b_years[4], 1, 0),
                            pars$term=="as.factor(year)2015"~ifelse(pars$conf.low <b_years[5] & pars$conf.high>b_years[5], 1, 0))

#Count number within for each parameter
count_within_95 <- aggregate(within_95 ~ term+model, pars, FUN=sum)
count_within_95$proportion <- count_within_95$within_95/length(fits)

##Bias
#Average error
mean_error <- aggregate(error ~ term+model, pars, FUN=mean)

##Combine into one table
pars_performance <- avg
pars_performance <- left_join(pars_performance, sd_avg, by=c("term", "model"))
pars_performance <- left_join(pars_performance, mean_error, by=c("term", "model"))
pars_performance <- left_join(pars_performance, error2, by=c("term", "model"))
pars_performance <- left_join(pars_performance, count_within_95, by=c("term", "model"))

#Save table
save(pars_performance, file="pars_performance.R")
write.csv(pars_performance, file="pars_performance.csv")

#Save pars

###Plots
#Distribution of mean error (to show bias)
#Boxplot
ggplot(pars, aes(y=error, x=term))+geom_boxplot(aes(fill=model, group=model))+facet_wrap("term", scales="free")+geom_hline(yintercept=0, linetype="dashed", linewidth=1.2)+theme(legend.position="none",  strip.text = element_blank())+
  facetted_pos_scales(
  y = list(NULL, ylim(-0.2,0.2), ylim(-0.15, 0.15), NULL, NULL, NULL, NULL, NULL, ylim(-0.4, 0.4), ylim(-0.15,0.15), ylim(-0.4, 0.4), ylim(-0.6, 0.6), ylim(-0.5, 0.5), ylim(-0.2, 0.2), ylim(-0.2, 0.2), NULL)
)
#Density plot w/ normal distribution curve
ggplot(pars, aes(x=error))+geom_density(aes(color=model, group=model))+facet_wrap("term", scales="free")+stat_function(fun = dnorm, args=list(mean = 0, sd = sd(pars$error)))+theme(legend.position="none")+xlim(-0.5, 0.5)

#Plot of number of simulations where true value falls within 95% confidence interval for each parameter for each value
ggplot(pars_performance, aes(x=model, y=proportion))+geom_line(aes(group=term, color=term))+theme(legend.position="none")+
  geom_text_repel(
    data = subset(pars_performance, pars_performance$model=="Data Rich\n(Correct\nprior)"),
    aes(label = term, color=term),
    size = 4,
    hjust=0.)