#******************************************************************
# Variable selection by Genetic Algorithm
# For the Matakana data
# This uses the data read by matakana_get_reference.R
# J. Dash (Jonathan.dash@scionresearch.com)
#******************************************************************

# Get libraries
library(GA)
library(reshape2)


###******************************************************************
#Calculate the control for model performance
###******************************************************************
control<-data.frame(me_BA=mean(ref$BasalArea),
                    me_TH=mean(ref$TopHeight),
                    me_SPH=mean(ref$Stocking),
                    me_TRV=mean(ref$TotalRecoverableVolume),
                    me_P=mean(ref$P_Sawlog),
                    me_UP=mean(ref$UP_Sawlog),
                    me_PULP=mean(ref$Pulp))

melt_cont<-melt(control)


ref_cont<-data.frame(c_ba=sqrt(mean((ref$BasalArea - control$me_BA)^2))/mean(ref$BasalArea),
                     c_TH=sqrt(mean((ref$TopHeight - control$me_TH)^2))/mean(ref$TopHeight),
                     c_sph=sqrt(mean((ref$Stocking - control$me_SPH)^2))/mean(ref$Stocking),
                     c_trv=sqrt(mean((ref$TotalRecoverableVolume - control$me_TRV)^2))/mean(ref$TotalRecoverableVolume),
                     c_P=sqrt(mean((ref$P_Sawlog - control$me_P)^2))/mean(ref$P_Sawlog),
                     c_UP=sqrt(mean((ref$UP_Sawlog - control$me_UP)^2))/mean(ref$UP_Sawlog),
                     c_pulp=sqrt(mean((ref$Pulp - control$me_PULP)^2))/mean(ref$Pulp))

ref_cont_melt<-melt(ref_cont)

# ################################################################################
# #Variable selection
# ###############################################################################
# 
# #Data frame containing all predictors 
# x<-subset(ref, select=-c( BasalArea, TopHeight, Stocking,  TotalRecoverableVolume, P_Sawlog, UP_Sawlog, Pulp
#                          ))
# 
# #Specify Y
# y<-ref[c("BasalArea", "TopHeight", "Stocking", "TotalRecoverableVolume", "P_Sawlog", "UP_Sawlog", "Pulp")]
# #y<-ref[c("Stocking")]
# #y<-ref[c("TopHeight")]
# 
# #Specify responses
# responses<-c("BasalArea", "TopHeight", "Stocking", "TotalRecoverableVolume", "P_Sawlog", "UP_Sawlog", "Pulp")
# #responses<-c("Stocking")
# #responses<-c("TopHeight")
# 
# #Imputation with all predictors
# set.seed(99)
# all<-yai(x=x, y=y, method="randomForest", k=1)
# yaiVarImp(all, plot=TRUE, nTop=30)
# 
# #write.csv(final, "final.csv")
# #write.csv(ref, "ref.csv")
# 
# #Summarise outputs
# grmsd(all, vars=responses)
# plot(all, vars=responses)
# ee<-rmsd(all, vars=responses, scale=F)
# 
# ###################################
# ##Variable selection with genetic algorithm
# ##################################
# fitness <- function(string) {
#   inc <- which(string == 1)
#   X <-  x[,inc]
#   mod <- yai(x=X, y=y, method="randomForest", k=1)
#   class(mod) <- "yai"
#   -grmsd(mod, vars=responses)
# }
# 
# #Select predictors using GA
# GA <- ga("binary", fitness = fitness, nBits = ncol(x),
#          names = colnames(x), monitor = plot, parallel=TRUE)
# 
# #plot(GA)
# summary(GA)
# 
# #Fit model with best subset of all predictors
# #set.seed(99)
# ga_imp<-yai(x=x[,GA@solution == 1], y=y, method="randomForest", k=1, ntree=10000)
# x2<-x[,GA@solution == 1]
# names(x2)
# ncol(x2)
# yaiVarImp(ga_imp, plot=TRUE, nTop=30)
# 
# #Summarise model outputs with GA
# grmsd(ga_imp, vars=responses)
# plot(ga_imp, vars=responses)
# dd<-rmsd(ga_imp, vars=responses, scale=F)
# 
# #Make a table to plot improvement for report
# #dd$RMSD_pct<-dd$rmsd
# dd<-cbind(dd, melt_cont[,2]) #Bring in the means
# colnames(dd)<-c("rmsd", "mean") #Rename to avoid confusion
# dd$rmsd_pct<-dd$rmsd/dd$mean #CAlculate RMSD as a percentage
# dd<-cbind(dd, ref_cont_melt[,2]) #Bring in the control RMSD as percentages
# colnames(dd) [4] <- "control_rmsd_pct" #Freshen the names
# dd$ratio<- dd$control_rmsd_pct/dd$rmsd_pct #Calcualte the ratios
# dd<-cbind(dd,ee)
# colnames(dd) [6] <- "all_pred__rmsd"
# dd$all_pred__rmsd_pct<-dd$all_pred__rmsd/dd$mean
# 
# 
# 
# #Further variable selection
# ###The plan is to use the outputs of the above as the inputs to the next GA
# fitness <- function(string) {
#   inc <- which(string == 1)
#   X <- x2[,inc]
#   mod <- yai(x=X, y=y, method="randomForest", k=1)
#   class(mod) <- "yai"
#   -grmsd(mod, vars=responses)
# }
# 
# #Select predictors using GA
# GA2 <- ga("binary", fitness = fitness, nBits = ncol(x2),
#           names = colnames(x2), monitor = plot, parallel=TRUE)
# summary(GA2)
# 
# x3<-x[,GA2@solution == 1]
# colnames(x3)
# ncol(x3)
# 
# x4<-x[,c("qav", "p30", "Cur_silvi_", "max.Ht.cr.base", "p70", "logcov", "ag.p99", "p50", "ag.p10")]
# 
# #Fit model with best subset of all predictors
# set.seed(999)
# ga_imp2<-yai(x=x4, y=y, method="randomForest", k=1, ntree=10000)
# grmsd(ga_imp2, vars=responses)
# plot(ga_imp2, vars=responses)
# rmsd(ga_imp2, vars=responses, scale=T)
# 
# yaiVarImp(ga_imp2, plot=T, nTop=30)
# 
# 
# head(ga_imp$neiDstRefs)
# ga_imputed<-impute(ga_imp, observed=T)
# hist(ga_imp$neiDstRefs, type="count")