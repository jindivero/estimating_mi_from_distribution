library(DHARMa)
library(viridis)
library(dplyr)
library(tidyr)
library(ggplot2)

#load fits
load("~/Dropbox/Mac/Documents/GitHub/thresholds_mi_distribution/thresholds_mi_distribution/current/output/model_fits.Rdata")

#load parameter estimates
pars_wide <- readRDS("~/Dropbox/Mac/Documents/GitHub/thresholds_mi_distribution/thresholds_mi_distribution/current/output/parameter_wide.rds")

##Spatial residual plot as percentile of Tweedie distribution

###DHARMa version
#Simulation-based residuals
simulate_fits <- function(x){
  s_res <- simulate(x, nsim = 500)
  return(s_res)
}

sims1 <- lapply(fits, simulate_fits)
sims2 <- lapply(fits_b, simulate_fits)
sims3 <- lapply(fits_c, simulate_fits)

#Dharma residuals
calculate_dharma <- function(sims, data_observed, data_predicted, column_preds, column_obs){
  observed <- as.numeric(data_observed[ , column_obs])
  predicted <- as.numeric(data_predicted[, column_preds])
  dharma <- DHARMa::createDHARMa(
    simulatedResponse = sims,
    observedResponse = observed,
    fittedPredictedResponse = predicted
  )
  return(dharma)
}

dharma1 <- mapply(FUN=calculate_dharma, sims1, data_sims_d1, preds1, "est_non_rf","observed", SIMPLIFY=F)
dharma2 <- mapply(FUN=calculate_dharma, sims2, data_sims_d1, preds2, "est_non_rf","observed", SIMPLIFY=F)
dharma3 <- mapply(FUN=calculate_dharma, sims3, data_sims_d1, preds3, "est_non_rf","observed", SIMPLIFY=F)

##Plot just one
dharma_test <- dharma1[[1]]
preds1_test <- preds1[[1]]
preds1_test$resids <- dharma_test$scaledResiduals
ggplot(preds1_test, aes(x=lon, y=lat))+geom_point(aes(color=resids), size=0.2)+scale_color_viridis()+facet_wrap("year")

plot(dharma_test)
DHARMa::testResiduals(dharma_test)




###Get quantiles from percentages using qtweedie--not a great method
#Function
calculate_quant <- function(data, column_probs, column_mus, ps, phis){
  probs <- as.numeric(data[ , column_probs])
  predicted <- as.numeric(data[, column_mus])
  quant <- qtweedie(p=probs, mu=predicted, power=ps, phi=phis)
  predicted2 <- cbind(data, quant)
  return(predicted2)
}

#Apply to each
q1 <- mapply(FUN=calculate_quant, nll1, "nll", "est", ps1, phis1, SIMPLIFY=F)
q2 <- mapply(FUN=calculate_quant, nll2, "nll", "est", ps2, phis2, SIMPLIFY=F)
q3 <- mapply(FUN=calculate_quant, nll3, "nll", "est", ps3, phis3, SIMPLIFY=F)

#Plot example of one
q_test <- q1[[1]]
ggplot(q_test, aes(x=lon, y=lat))+geom_point(aes(color=quant), size=0.2)+scale_color_viridis()+facet_wrap("year")
