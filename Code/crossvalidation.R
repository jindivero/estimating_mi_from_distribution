###Cross-fold validation
library(cvAUC)
library(dplyr)
library(tidyr)
library(rlang)
library(purrr)
library(mgcv)
library(Metrics)
library(viridis)
install.packages("cvAUC")
library(ggplot2)
install.packages("tweedie")
library(tweedie)
install.packages("poistweedie")
library(poistweedie)
library(DHARMa)
library(KernSmooth)
library(graphics)
library(MASS)

##Set ggplot themes
theme_set(theme_bw(base_size = 15))
theme_update(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#load simulated data
data_sims <- readRDS("~/Dropbox/Mac/Documents/GitHub/thresholds_mi_distribution/thresholds_mi_distribution/current/output/data_sims.rds")

#Remove mesh from data list
data_sims_d <- flatten(data_sims)
data_sims_d1 <- keep(.x=data_sims_d, .p=is.data.frame)

#load cross validation data
data_sims_CV <- readRDS("~/Dropbox/Mac/Documents/GitHub/thresholds_mi_distribution/thresholds_mi_distribution/current/output/data_sims_CV.rds")

##Remove mesh from data list
data_sims_CV_d <- flatten(data_sims_CV)
data_sims_CV_d1 <- keep(.x=data_sims_CV_d, .p=is.data.frame)

#Load predictions
load("~/Dropbox/Mac/Documents/GitHub/thresholds_mi_distribution/thresholds_mi_distribution/current/output/preds.Rdata")

#load parameter estimates
pars_wide <- readRDS("~/Dropbox/Mac/Documents/GitHub/thresholds_mi_distribution/thresholds_mi_distribution/current/output/parameter_wide.rds")

#load fits
load("~/Dropbox/Mac/Documents/GitHub/thresholds_mi_distribution/thresholds_mi_distribution/current/output/model_fits.Rdata")

#Extract tweedie parameters (phi's and p's)
phis1 <- as.matrix(pars_wide[c(1:100), 16])
phis2<- as.matrix(pars_wide[c(101:200), 16])
phis3 <- as.matrix(pars_wide[c(201:300), 16])
ps1<- as.matrix(pars_wide[c(1:100), 18])
ps2 <- as.matrix(pars_wide[c(101:200), 18])
ps3<- as.matrix(pars_wide[c(201:300), 18])

#Negative log probability of observations given model fit

calculate_nll <- function(data_observed, data_predicted, column_obs, column_preds, ps, phis){
  observed <- as.numeric(data_observed[ , column_obs])
  predicted <- as.numeric(data_predicted[, column_preds])
  nll <- dtweedie(y=observed, mu=predicted,power=ps, phi=phis)
  predicted2 <- cbind(data_predicted, nll)
  return(predicted2)
}

#Apply to each
nll1 <- mapply(FUN=calculate_nll, data_sims_CV_d1, preds_CV1, "observed", "est", ps1, phis1, SIMPLIFY=F)
nll2 <- mapply(FUN=calculate_nll, data_sims_CV_d1, preds_CV2, "observed", "est", ps2, phis2, SIMPLIFY=F)
nll3 <- mapply(FUN=calculate_nll, data_sims_CV_d1, preds_CV3, "observed", "est", ps3, phis3, SIMPLIFY=F)

#Sum overall
sum_nll <- function(nll, column_nll){
  nlls <- as.numeric(nll[ , column_nll])
  sum_nll <- sum(nlls)
  return(sum_nll)
}
nll_sum1 <- mapply(FUN=sum_nll, nll1, "nll", SIMPLIFY=F)
nll_sum2 <- mapply(FUN=sum_nll, nll2, "nll", SIMPLIFY=F)
nll_sum3 <- mapply(FUN=sum_nll, nll3, "nll", SIMPLIFY=F)

#Sum in areas/times where below simulated true threshold (=s95)
#
delta <- 2
x50 <- 2
s95 <- x50+delta


#Subset to below s50
nll1_thresh <- lapply(nll1, subset, mi_0.3<2)
nll2_thresh <- lapply(nll2, subset, mi_0.3<2)
nll3_thresh <- lapply(nll3, subset, mi_0.3<2)

#Subset to where observations have no fish
#nll1_thresh <- lapply(nll1, subset, pred2==0)
#nll2_thresh <- lapply(nll2, subset, pred2==0)
#nll3_thresh <- lapply(nll3, subset, pred2==0)

nll_sum_thresh1 <- mapply(FUN=sum_nll, nll1_thresh, "nll", SIMPLIFY=F)
nll_sum_thresh2 <- mapply(FUN=sum_nll, nll2_thresh, "nll", SIMPLIFY=F)
nll_sum_thresh3 <- mapply(FUN=sum_nll, nll3_thresh, "nll", SIMPLIFY=F)

#Subset to above s95
nll1_95 <- lapply(nll1, subset, mi_0.3>4)
nll2_95 <- lapply(nll2, subset, mi_0.3>4)
nll3_95 <- lapply(nll3, subset, mi_0.3>4)

nll1_sum_95 <- mapply(FUN=sum_nll, nll1_95, "nll", SIMPLIFY=F)
nll2_sum_95 <- mapply(FUN=sum_nll, nll2_95, "nll", SIMPLIFY=F)
nll3_sum_95 <- mapply(FUN=sum_nll, nll3_95, "nll", SIMPLIFY=F)


##Combine
#Overall
nll_all <- as.data.frame(unlist(nll_sum1))
nll_all$model2 <- as.vector(unlist(nll_sum2))
nll_all$model3 <- as.vector(unlist(nll_sum3))
colnames(nll_all) <- c("No Prior", "Typical Case", "Unusual Case")

#Above s95
nll_95 <- as.data.frame(unlist(nll1_sum_95))
nll_95$model2 <- as.vector(unlist(nll2_sum_95))
nll_95$model3 <- as.vector(unlist(nll3_sum_95))
colnames(nll_95) <- c("No Prior", "Typical Case", "Unusual Case")

#Below s95
nll_all_thresh <- as.data.frame(unlist(nll_sum_thresh1))
nll_all_thresh$model2 <- as.vector(unlist(nll_sum_thresh2))
nll_all_thresh$model3 <- as.vector(unlist(nll_sum_thresh3))
colnames(nll_all_thresh) <- c("No Prior", "Typical Case", "Unusual Case")

#Identify minimum in each row
nll_all$min <- as.numeric(apply(nll_all, 1, FUN = min))
nll_all_thresh$min <- as.numeric(apply(nll_all_thresh, 1, FUN = min))
nll_95$min <- as.numeric(apply(nll_95, 1, FUN = min))

#Calculate difference
nll_all$No_Prior_diff <- nll_all$`No Prior`-nll_all$min
nll_all$Typical_Case_diff <- nll_all$`Typical Case`-nll_all$min
nll_all$Unusual_Case_diff <- nll_all$`Unusual Case`-nll_all$min

nll_all_thresh$No_Prior_diff <- nll_all_thresh$`No Prior`-nll_all_thresh$min
nll_all_thresh$Typical_Case_diff <- nll_all_thresh$`Typical Case`-nll_all_thresh$min
nll_all_thresh$Unusual_Case_diff <- nll_all_thresh$`Unusual Case`-nll_all_thresh$min

nll_95$No_Prior_diff <- nll_95$`No Prior`-nll_95$min
nll_95$Typical_Case_diff <- nll_95$`Typical Case`-nll_95$min
nll_95$Unusual_Case_diff <- nll_95$`Unusual Case`-nll_95$min

#Add all sim
nll_all_thresh$sim <- 1:100
nll_all$sim <- 1:100
nll_95$sim <- 1:100

#Combine for plotting
#Add column for type
#Rename columns
nll_all$Type <-"Overall"
nll_95$Type <- "Above_s95"
nll_all_thresh$Type <-"Below_s50"

#Remove min
nll_all$min <- NULL
nll_all_thresh$min <- NULL
nll_95$min <- NULL

#Pivot long
nll_all <- pivot_longer(nll_all, c(1:6), names_to="Model")
nll_all_thresh <- pivot_longer(nll_all_thresh, c(1:6), names_to="Model")
nll_95 <- pivot_longer(nll_95, c(1:6), names_to="Model")

#Relabel
colnames(nll_all)[4] <- "Overall"
colnames(nll_all_thresh)[4] <- "Below_s50"
colnames(nll_95)[4] <- "Above_s95"

nll_all$Type <- NULL
nll_all_thresh$Type <- NULL
nll_95$Type <- NULL

#Combine
nll_combined <- merge(nll_95, nll_all_thresh, by=c("sim", "Model"))

#Separate difference
nll_combined_diff <- subset(nll_combined, Model=="No_Prior_diff"|Model=="Typical_Case_diff"|Model=="Unusual_Case_diff")

#Separate raw values
nll_combined_raw <- subset(nll_combined, Model=="No Prior"|Model=="Typical Case"|Model=="Unusual Case")

####Plot density Fried egg (aka a haunted house plot)
ggplot(nll_combined_diff, aes(x=Above_s95, y=Below_s50))+
  stat_density_2d(h = c(ifelse(bandwidth.nrd(nll_combined_diff$Above_s95) == 0, 0.1, bandwidth.nrd(nll_combined_diff$Above_s95)),ifelse(bandwidth.nrd(nll_combined_diff$Below_s50) == 0, 0.1, bandwidth.nrd(nll_combined_diff$Below_s50))),
                  geom = "raster",aes(fill = after_stat(density)),contour = FALSE) + 
  scale_fill_viridis_c()+geom_point(color="white", size=3.5, shape="*")+
  facet_grid("Model")

ggplot(nll_combined_raw, aes(x=Above_s95, y=Below_s50))+
  stat_density_2d(geom = "raster",aes(fill = after_stat(density)),contour = FALSE)+
  scale_fill_viridis_c()+geom_point(color="white", size=3.5, shape="*")+
  facet_grid("Model")+
  xlab("Above s95")+
  ylab("Below s50")


#ggplot




#Other plot options
ggplot(nll_combined_diff, aes(x=Model, y=Below_Threshold))+geom_boxplot()
ggplot(nll_combined_diff, aes(x=Model, y=Overall))+geom_boxplot()

###Tim's code

# first create vectors of the two things to be plotted - here they were "estimated.beta" and "estimated.sd.r".  In your case it will be the two performance measures.
overall1 <- nll_combined_diff$Overall
below1 <- nll_combined_diff$Below_Threshold

kernel.dens<-bkde2D(cbind(overall1,below1),bandwidth=c(0.001, 0.001),range.x=list(c(-0.0,0.08),c(0.0,0.08)),truncate=FALSE,gridsize=c(2000,2000))
paletteable.egg<-colorRampPalette(c("#BFEFFF","white","white", "yellow","#FFC125"))
xlims<-c(-0.0,0.1)
ylims<-c(0.0,0.1)

filled.contour(kernel.dens$x1, kernel.dens$x2, kernel.dens$fhat, nlevels=15, color.palette =paletteable.egg,
                   xlab="", ylab="", xlim=xlims, ylim=ylims, cex.lab=2)
x.lab.text<-"Overall"
y.lab.text<-"Below Threshold"

box()
axis(1,cex.axis=1.5)
axis(2,cex.axis=1.5)
mtext(side=1,x.lab.text,line=3,cex=2)
mtext(side=2,y.lab.text,line=3,cex=2)

par(new=T,xpd=NA)
points(overall,below,pch=19,col="black")


#####OLDER VERSIONS
###Calculate root mean squared error
##Create function to calculate RMSE 
calculate_rmse <- function(data_observed, data_predicted, column_obs, column_preds){
  observed <- as.numeric(data_observed[ , column_obs])
  predicted <- as.numeric(data_predicted[, column_preds])
  rmse1 <- rmse(actual=observed, predicted=predicted)
  return(rmse1)
}

##Apply for cross-validation
rmse_CV1<- mapply(FUN=calculate_rmse, data_sims_CV_d1, preds_CV1, "observed","pred2", SIMPLIFY=F)
rmse_CV2<- mapply(FUN=calculate_rmse, data_sims_CV_d1, preds_CV2, "observed","pred2", SIMPLIFY=F)
rmse_CV3<- mapply(FUN=calculate_rmse, data_sims_CV_d1, preds_CV3, "observed","pred2", SIMPLIFY=F)

##Apply for model residuals
rmse1<- mapply(FUN=calculate_rmse, data_sims_d1, preds1, "observed","pred2", SIMPLIFY=F)
rmse2<- mapply(FUN=calculate_rmse, data_sims_d1, preds2, "observed","pred2", SIMPLIFY=F)
rmse3<- mapply(FUN=calculate_rmse, data_sims_d1, preds3, "observed","pred2", SIMPLIFY=F)

###Plot model residuals for plotting
##Combine into one
rmse_all <- as.data.frame(unlist(rmse1))
rmse_all$model2 <- as.vector(unlist(rmse2))
rmse_all$model3 <- as.vector(unlist(rmse3))
colnames(rmse_all) <- c("Data Limited", "Data Rich", "Unusual Case")
rmse_all$sim <- 1:100
#Pivot long
rmse_all <- pivot_longer(rmse_all, 1:3, names_to="Model")

##Plot
#Boxplot
ggplot(rmse_all, aes(x=Model, y=value, group=Model, color=Model))+geom_boxplot()+theme(legend.position="none")

#Line plot
ggplot(rmse_all, aes(x=sim, y=value, group=Model, color=Model))+geom_line()

###For CVs
rmse_all_CV <- as.data.frame(unlist(rmse_CV1))
rmse_all_CV$model2 <- as.vector(unlist(rmse_CV2))
rmse_all_CV$model3 <- as.vector(unlist(rmse_CV3))
colnames(rmse_all_CV) <- c("Data Limited", "Data Rich", "Unusual Case")
rmse_all_CV$sim <- 1:100
#Pivot long
rmse_all_CV <- pivot_longer(rmse_all_CV, 1:3, names_to="Model")

##Plot
#Boxplot
ggplot(rmse_all_CV, aes(x=Model, y=value, group=Model, color=Model))+geom_boxplot()+theme(legend.position="none")

#Line plot
ggplot(rmse_all, aes(x=sim, y=value, group=Model, color=Model))+geom_line()+theme(legend.position="none")+ylab("Root Mean Square Error")

###Difference in RMSE, instead of raw value
##Recreate wide dataset
rmse_all_CV <- as.data.frame(unlist(rmse_CV1))
rmse_all_CV$model2 <- as.vector(unlist(rmse_CV2))
rmse_all_CV$model3 <- as.vector(unlist(rmse_CV3))
colnames(rmse_all_CV) <- c("Data Limited", "Data Rich", "Unusual Case")
#Identify minimum RMSE in a row
rmse_all_CV$min <- apply(rmse_all_CV, 1, FUN = min)
#Calculate difference
rmse_all_CV$Data_Limited <- rmse_all_CV$`Data Limited`-rmse_all_CV$min
rmse_all_CV$Data_Rich <- rmse_all_CV$`Data Rich`-rmse_all_CV$min
rmse_all_CV$Unusual_Case <- rmse_all_CV$`Unusual Case`-rmse_all_CV$min
#add sim column
rmse_all_CV$sim <- 1:100
#Pivot longer again
rmse_all_CV <- pivot_longer(rmse_all_CV[c(5:8)], 1:3, names_to="Model")
#Plot
ggplot(rmse_all_CV, aes(x=Model, y=value, group=Model, color=Model))+geom_boxplot()+theme(legend.position="none")+ylab("Delta RMSE")

###Map of prediction residuals
##Just one
dat_test <- data_sims_d1[[1]]
pred_test <- preds1[[1]]
dat_test$residual <- pred_test$pred2-dat_test$observed
ggplot(dat_test, aes(x=lon, y=lat))+geom_point(aes(color=residual), size=0.4)+scale_color_viridis()+facet_wrap("year")

####Comparing to presence-absence 
###Create a presence-absence column
calculate_auc <- function(data_observed,data_predicted, col_obs, col_preds){
  predicted <- as.numeric(ifelse(data_predicted[,col_preds] >0,1,0))
  observed <- as.numeric(ifelse(data_observed[,col_obs] >0,1,0))
  auc <- auc(actual=observed, predicted=predicted)
  return(auc)
}

##
auc1<- mapply(FUN=calculate_auc, data_sims_CV_d1, preds1, "observed","pred2", SIMPLIFY=F)
auc2<- mapply(FUN=calculate_auc, data_sims_CV_d1, preds2, "observed","pred2", SIMPLIFY=F)
auc3<- mapply(FUN=calculate_auc, data_sims_CV_d1, preds3, "observed","pred2", SIMPLIFY=F)

#Plot
auc_all_CV <- as.data.frame(unlist(auc1))
auc_all_CV$model2 <- as.vector(unlist(auc2))
auc_all_CV$model3 <- as.vector(unlist(auc3))
colnames(auc_all_CV) <- c("Data Limited", "Data Rich", "Unusual Case")
auc_all_CV$sim <- 1:100
#Pivot long
auc_all_CV <- pivot_longer(auc_all_CV, 1:3, names_to="Model")

##Plot
#Boxplot
ggplot(auc_all_CV, aes(x=Model, y=value, group=Model, color=Model))+geom_boxplot()+theme(legend.position="none")+ylab("Area Under the Receiver Operating Curve")





#Threshold testing
test$mi_effect1 <- 4 * (1 / (1 + exp(-log(19) * (test$mi_0.3 - 2) / 2)) - 1)
ggplot(test, aes(x=mi_0.3, y=mi_effect1))+geom_point()
ggplot(test,aes(x=mi_0.3))+geom_density(color="purple")+geom_density(test,mapping=aes(x=invtemp), color="red")+geom_density(test,mapping=aes(x=po2), color="blue")+xlab("MI=purple, pO2=blue, invtemp=red")
ggplot(test,aes(x=invtemp))+geom_density()
ggplot(test,aes(x=po2))+geom_density()
