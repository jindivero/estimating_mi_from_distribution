
###Extract convergence diagnostics
##Function
extract_convergence <- function(x){
  if(!is.character(x)){
    convergence <- as.data.frame(x[["model"]][["convergence"]])
    colnames(convergence) <- "convergence"
    convergence$iterations <- x$model$iterations
    convergence$message <- x$model$message
    convergence$evaluations_function <- x$model$evaluations[1]
    convergence$evaluations_gradient <- x$model$evaluations[2]
  }
  if(is.character(x)){
    convergence <- as.data.frame(matrix(ncol=5, nrow=1, NA))
    colnames(convergence) <- c("convergence", "iterations", "message", "evaluations_function", "evaluations_gradient")
    convergence[,3] <- paste(x[[1]])
  }
  cons <- bind_rows(convergence)
  return(cons)
}
##Apply to each model
convergence1 <- lapply(fits, extract_convergence)
convergence2 <- lapply(fits2, extract_convergence)
convergence3 <- lapply(fits3, extract_convergence)
convergence4 <- lapply(fits4, extract_convergence)
convergence5 <- lapply(fits5, extract_convergence)

##Bind together
convergence1 <- bind_rows(convergence1)
convergence2 <- bind_rows(convergence2)
convergence3 <- bind_rows(convergence3)
convergence4 <- bind_rows(convergence4)
convergence5 <- bind_rows(convergence5)

##Add column for data simulation and model
convergence1$sim <- 1:100
convergence1$model <- model_names[1]
convergence2$sim <- 1:100
convergence2$model <- model_names[2]
convergence3$sim <- 1:100
convergence3$model <- model_names[3]
convergence4$sim <- 1:100
convergence4$model <- names(fits_all[4])
convergence5$sim <- 1:100
convergence5$model <- names(fits_all[5])

#Bind all together
convergence <-bind_rows(convergence1, convergence2, convergence3, convergence4) #convergence5)
#Truncate message
convergence$message <- substr(convergence$message, 1, 100)
##Summarize number of each type of convergence for each model
cons <- convergence %>% group_by(message, model) %>% summarise(count=n()) 
##Flip to wide
cons <-pivot_wider(cons, names_from=model, values_from=count)
colnames(cons)[1] <-"metric"

###Models that ran without error
fitsa <- purrr::discard(fits, is.character)
fits2a <- purrr::discard(fits2, is.character)
fits3a <- purrr::discard(fits3, is.character)
fits4a <- purrr::discard(fits4, is.character)
fits5a <- purrr::discard(fits5, is.character)

###Make dataframe counting number and add
fits_count <- as.data.frame(matrix(nrow=1, ncol=6))
colnames(fits_count) <- c("metric", model_names)
fits_count[1,1] <- "Completed"
fits_count[1,2] <- length(fitsa)
fits_count[1,3] <- length(fits2a)
fits_count[1,4] <- length(fits3a)
fits_count[1,5] <- length(fits4a)
fits_count[1,6] <- length(fits5a)

###Combine convergence and completed models
diagnostics <- bind_rows(cons, fits_count)

###Positive definite hessian
##Function to extract for each simulation
extract_pdHess <- function(x){
  if(!is.character(x)){
    pdh <- as.data.frame(x$pos_def_hessian)
    return(pdh)
  }
}
##Apply to each model
pdHess1 <- lapply(fits, extract_pdHess)
pdHess2 <- lapply(fits2, extract_pdHess)
pdHess3 <- lapply(fits3, extract_pdHess)
pdHess4 <- lapply(fits4, extract_pdHess)
pdHess5 <- lapply(fits5, extract_pdHess)

##Function to get sum of number of simulations with non-positive definite hessian in each model
extract_pdHess2 <- function(pdHess){
  pdHess <- bind_rows(pdHess)
  pdHess$sim <- as.numeric(row.names(pdHess))
  pdHess$sim <- as.character(pdHess$sim)
  pdHess$neg_def_hessian <- ifelse(pdHess$"x$pos_def_hessian"=="TRUE", 1,0)
  number_negdHess <- print(sum(pdHess$neg_def_hessian))
  return(number_negdHess)
}

##Make dataframe to store
hess <- as.data.frame(matrix(nrow=1, ncol=6))
colnames(hess) <- c("metric", model_names)
hess[1,1] <- "Positive definite hessian matrix"
hess[1,2] <- extract_pdHess2(pdHess1)
hess[1,3] <- extract_pdHess2(pdHess2)
hess[1,4] <- extract_pdHess2(pdHess3)
hess[1,5] <- extract_pdHess2(pdHess4)
hess[1,6] <- extract_pdHess2(pdHess5)


###Maximum final gradient on each parameter
##Function to exactract gradients from each model
extract_grad <- function(x){
  if(!is.character(x)){
    grad <- as.data.frame(x$gradients)
    return(grad)
  }
}
##Apply to all models
grad1<- lapply(fits, extract_grad)
grad2<- lapply(fits2, extract_grad)
grad3<- lapply(fits3, extract_grad)
grad4<- lapply(fits4, extract_grad)
grad5<- lapply(fits5, extract_grad)

clean_grad <- function(grad, fits){
  if(!is.null(grad)){
    grad <- bind_cols(grad)
  }
  grad <- as.data.frame(t(grad))
  names <- names(fits[[1]]$sd_report$par.fixed)
  colnames(grad) <- pars_names
  grad$sim <- 1:nrow(grad)
  grad <- pivot_longer(grad,cols=1:11, names_to="term", values_to="gradient")
  grad$big_grad <- ifelse(grad$gradient<0.001, 1,0)
  large_grad <- aggregate(big_grad ~ term, grad, FUN=sum)
  #large_grad$proportion <- large_grad$big_grad/length(fits)
  return(large_grad)
}

##Run for each model
grads1<- clean_grad(grad1, fits)
grads2<- clean_grad(grad2, fits2)
grads3<- clean_grad(grad3, fits3)
grads4<- clean_grad(grad4, fits4)
grads5<- clean_grad(grad5, fits5)

##Bind together
grads1$model2 <-grads2$big_grad
grads1$model3 <-grads3$big_grad
grads1$model4 <-grads4$big_grad
grads1$model5 <-grads5$big_grad
colnames(grads1) <-c("metric", model_names)

####Combine gradient into final diagnostics
diagnostics <- bind_rows(diagnostics, grads1)

###Plots
##Convergence
##Flip to long for plotting
cons_long <- pivot_longer(cons, cols=2:6, names_to="model")
cons_long$model <- str_wrap(cons_long$model, width=10)
cons_long$metric <- str_wrap(cons_long$metric, width=40)
ggplot(cons_long, aes(x=model, y=value))+geom_col(aes(group=metric, fill=metric))







###~~~~~~~~~~~~Looking closer at depths and spatial effects

#flatten predictions into one dataframe
preds_mi2 <-bind_rows(preds_mi)
preds_mi2$rf <- preds_mi2$est-preds_mi2$est_non_rf

ggplot(preds_mi2, aes(x=omega_s, y=log_depth_sc2))+geom_point()

#Look at correlation of parameter estimates
ggplot(pars_wide, aes(x=log_depth_sc2, y=range))+geom_point()+geom_hline(yintercept=range)+geom_vline(xintercept=beta2)
ggplot(pars_wide, aes(x=depth_sc2, y=sigma_O))+geom_point()+geom_hline(yintercept=sigma_O)+geom_vline(xintercept=beta2)

#Look at spatial random effects across space
ggplot(preds_test, aes(lat, lon, col = omega_s)) +
  scale_colour_gradient2() +
  geom_point(size=0.5)

ggplot(preds_test, aes(x=omega_s, y=depth_scaled2))+geom_point()

