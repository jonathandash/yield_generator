#*****************************************************************************
# Variable selection using the top 5 approach for Matakana.
# This should be made flexible eventually
# J. Dash (jonathan.dash@scionresearch.com) 
# January 2016
################################################################################

library(yaImpute)

#Data frame containing all predictors 
x<-subset(ref, select=-c( BasalArea, TopHeight, Stocking,  TotalRecoverableVolume, P_Sawlog, UP_Sawlog, Pulp, Source
))

#Specify Y
y_PS<-ref[c("P_Sawlog")]

#Imputation with all predictors
P.imp<-yai(x=x, y=y_PS, method="randomForest", k=5, ntree=100000)

png(paste(output.dir, "Pruned_VarbImportance.png", sep=""), w=20, h=20, units="cm", res=300)
yaiVarImp(P.imp, plot=TRUE, nTop=40)
dev.off()

png(paste(output.dir, "Pruned_Imputation.png", sep=""), w=20, h=20, units="cm", res=300)
plot(P.imp, vars= "P_Sawlog")
dev.off()

rmsd(P.imp, vars="P_Sawlog")

#plot(max.Ht.cr.base ~ P_Sawlog, data=ref)
#boxplot(P_Sawlog ~ Cur_silvi_, data=ref)

ref.crptyp<-ref %>% group_by(Cur_silvi_) %>% summarise(mean_P_Sawlog = mean(P_Sawlog), 
                                                       min_p_Sawlog=min(P_Sawlog),
                                                       max_P_Sawlog=max(P_Sawlog),
                                                       count = length(P_Sawlog))

# png(paste(output.dir, "P_Sawlog.png", sep=""), w=20, h=20, units="cm", res=300)
# ggplot (ref, aes(Age, P_Sawlog)) +geom_point(aes(colour=Cur_silvi_))
# dev.off()

#Optimise k
resultsPS<-data.frame(k=numeric(), rmse=numeric())
for (i in 1:20)
{
  set.seed(99)
  t<-yai(x=x, y=y_PS, method="randomForest", k=i)
  timp<-impute(t, method="dstWeighted", k=i, observed=T)
  timp$err<-timp$P_Sawlog - timp$P_Sawlog.o
  rmse<-sqrt(mean(timp$err^2))
  resultsPS<-rbind(resultsPS, c(i, rmse))
}

colnames(resultsPS)<- c('k', 'rmse')
resultsPS<-resultsPS[with(resultsPS, order(rmse)), ] #Order by RMSE
print(paste("Optimal value for k is ", resultsPS [1,1], sep=""))

#Model with optimal k
#PS.imp<-yai(x=x, y=y_PS, method="randomForest", k=resultsPS[1,1], ntree=100000)
#PS.imp.out<-yai(x=x, y=y_PS, method="P_Sawlog", observed=TRUE, method='dstWeighted')


#Specify Y
y_TRV<-ref[c("TotalRecoverableVolume")]

#Imputation with all predictors for TRV
TRV.imp<-yai(x=x, y=y_TRV, method="randomForest", k=5, ntree=100000)

png(paste(output.dir, "TRV_VarbImportance.png", sep=""), w=20, h=20, units="cm", res=300)
yaiVarImp(TRV.imp, plot=TRUE, nTop=40)
dev.off()

png(paste(output.dir, "TRV_imputation.png", sep=""), w=20, h=20, units="cm", res=300)
plot(TRV.imp, vars= "TotalRecoverableVolume")
dev.off()

#Optimise k
resultsTRV<-data.frame(k=numeric(), rmse=numeric())
for (i in 1:20)
{
  set.seed(99)
  t<-yai(x=x, y=y_TRV, method="randomForest", k=i)
  timp<-impute(t, method="dstWeighted", k=i, observed=T)
  timp$err<-timp$TotalRecoverableVolume - timp$TotalRecoverableVolume.o
  rmse<-sqrt(mean(timp$err^2))
  resultsTRV<-rbind(resultsTRV, c(i, rmse))
}

colnames(resultsTRV)<- c('k', 'rmse')
resultsTRV<-resultsTRV[with(resultsTRV, order(rmse)), ] #Order by RMSE
print(paste("Optimal value for k is ", resultsTRV [1,1], sep=""))

#Model with optimal k
#TRV.imp<-yai(x=x, y=y_TRV, method="randomForest", k=resultsTRV[1,1], ntree=100000)
#TRV.imp.out<-yai(x=x, y=y_TRV, method="TotalRecoverableVolume", observed=TRUE, method='dstWeighted')

###
#Stocking
###

#Specify Y
y_SPH<-ref[c("Stocking")]

#Imputation with all predictors for TRV
SPH.imp<-yai(x=x, y=y_SPH, method="randomForest", k=5, ntree=100000)

png(paste(output.dir, "SPH_VarbImportance.png", sep=""), w=20, h=20, units="cm", res=300)
yaiVarImp(SPH.imp, plot=TRUE, nTop=40)
dev.off()

png(paste(output.dir, "SPH_Imputation.png", sep=""), w=20, h=20, units="cm", res=300)
plot(SPH.imp, vars= "Stocking")
dev.off()

#Optimise k
resultsSPH<-data.frame(k=numeric(), rmse=numeric())
for (i in 1:20)
{
  set.seed(99)
  t<-yai(x=x, y=y_SPH, method="randomForest", k=i)
  timp<-impute(t, method="dstWeighted", k=i, observed=T)
  timp$err<-timp$Stocking - timp$Stocking.o
  rmse<-sqrt(mean(timp$err^2))
  resultsSPH<-rbind(resultsSPH, c(i, rmse))
}

colnames(resultsSPH)<- c('k', 'rmse')
resultsSPH<-resultsSPH[with(resultsSPH, order(rmse)), ] #Order by RMSE
print(paste("Optimal value for k is ", resultsSPH [1,1], sep=""))

#Model with optimal k
#SPH.imp<-yai(x=x, y=y_SPH, method="randomForest", k=resultsSPH[1,1], ntree=100000)
#SPH.imp.out<-yai(x=x, y=y_SPH, method="Stocking", observed=TRUE, method='dstWeighted')

###
#Basal Area
###

#Specify Y
y_BA<-ref[c("BasalArea")]

#Imputation with all predictors for TRV
BA.imp<-yai(x=x, y=y_BA, method="randomForest", k=5, ntree=100000)

png(paste(output.dir, "BA_VarbImportance.png", sep=""), w=20, h=20, units="cm", res=300)
yaiVarImp(BA.imp, plot=TRUE, nTop=40)
dev.off()

png(paste(output.dir, "BA_Imputation.png", sep=""), w=20, h=20, units="cm", res=300)
plot(BA.imp, vars= "BasalArea")
dev.off()

#Optimise k
resultsBA<-data.frame(k=numeric(), rmse=numeric())
for (i in 1:20)
{
  set.seed(99)
  t<-yai(x=x, y=y_BA, method="randomForest", k=i)
  timp<-impute(t, method="dstWeighted", k=i, observed=T)
  timp$err<-timp$BasalArea - timp$BasalArea.o
  rmse<-sqrt(mean(timp$err^2))
  resultsBA<-rbind(resultsBA, c(i, rmse))
}

colnames(resultsBA)<- c('k', 'rmse')
resultsBA<-resultsBA[with(resultsBA, order(rmse)), ] #Order by RMSE
print(paste("Optimal value for k is ", resultsBA [1,1], sep=""))

#Model with optimal k
#BA.imp<-yai(x=x, y=y_BA, method="randomForest", k=resultsBA[1,1], ntree=100000)
#BA.imp.out<-impute(BA.imp, vars="BasalArea", observed=TRUE, method='dstWeighted')

###
#UP_Sawlog
###

#Specify Y
y_UP<-ref[c("UP_Sawlog")]

#Imputation with all predictors for UP
UP.imp<-yai(x=x, y=y_UP, method="randomForest", k=5, ntree=100000)

png(paste(output.dir, "UP_VarbImportance.png", sep=""), w=20, h=20, units="cm", res=300)
yaiVarImp(UP.imp, plot=TRUE, nTop=40)
dev.off()

png(paste(output.dir, "UP_Imputation.png", sep=""), w=20, h=20, units="cm", res=300)
plot(UP.imp, vars= "UP_Sawlog")
dev.off()

#Optimise k
resultsUP<-data.frame(k=numeric(), rmse=numeric())
for (i in 1:20)
{
  set.seed(99)
  t<-yai(x=x, y=y_UP, method="randomForest", k=i)
  timp<-impute(t, method="dstWeighted", k=i, observed=T)
  timp$err<-timp$UP_Sawlog - timp$UP_Sawlog.o
  rmse<-sqrt(mean(timp$err^2))
  resultsUP<-rbind(resultsUP, c(i, rmse))
}
colnames(resultsUP)<- c('k', 'rmse')
resultsUP<-resultsUP[with(resultsUP, order(rmse)), ] #Order by RMSE
print(paste("Optimal value for k is ", resultsUP [1,1], sep=""))

#Model with optimal k
#UP.imp<-yai(x=x, y=y_UP, method="randomForest", k=resultsUP[1,1], ntree=100000)
#UP.imp.out<-impute(UP.imp, vars="UP_Sawlog", observed=TRUE, method='dstWeighted')

#******
##Pulp##
#*****

#Specify Y
y_Pulp<-ref[c("Pulp")]

#Imputation with all predictors for TRV
Pulp.imp<-yai(x=x, y=y_Pulp, method="randomForest", k=5, ntree=100000)

png(paste(output.dir, "Pulp_VarbImportance.png", sep=""), w=20, h=20, units="cm", res=300)
yaiVarImp(Pulp.imp, plot=TRUE, nTop=40)
dev.off()

png(paste(output.dir, "Pulp_Imputation.png", sep=""), w=20, h=20, units="cm", res=300)
plot(Pulp.imp, vars= "Pulp")
dev.off()
#Pulp.imp.out<-impute(Pulp.imp, vars="Pulp", observed=TRUE, method='dstWeighted')

#Optimise k
resultsPu<-data.frame(k=numeric(), rmse=numeric())
for (i in 1:20)
{
  set.seed(99)
  t<-yai(x=x, y=y_Pulp, method="randomForest", k=i)
  timp<-impute(t, method="dstWeighted", k=i, observed=T)
  timp$err<-timp$Pulp - timp$Pulp.o
  rmse<-sqrt(mean(timp$err^2))
  resultsPu<-rbind(resultsPu, c(i, rmse))
}
colnames(resultsPu)<- c('k', 'rmse')
resultsPu<-resultsPu[with(resultsPu, order(rmse)), ] #Order by RMSE
print(paste("Optimal value for k is ", resultsPu [1,1], sep=""))

#Model with optimal k
#Pulp.imp<-yai(x=x, y=y_Pulp, method="randomForest", k=resultsPu[1,1], ntree=100000)
#Pulp.imp.out<-impute(Pulp.imp, vars="Pulp", observed=TRUE, method='dstWeighted')

#***
##TopHeight##
#***

#Specify Y
y_TH<-ref[c("TopHeight")]

#Imputation with all predictors for TRV
TH.imp<-yai(x=x, y=y_TH, method="randomForest", k=5, ntree=100000)

png(paste(output.dir, "TH_VarbImportance.png", sep=""), w=20, h=20, units="cm", res=300)
yaiVarImp(TH.imp, plot=TRUE, nTop=40)
dev.off()

png(paste(output.dir, "TH_Imputation.png", sep=""), w=20, h=20, units="cm", res=300)
plot(TH.imp, vars= "TopHeight")
dev.off()

#Optimise k
results99<-data.frame(k=numeric(), rmse=numeric())
for (i in 1:20)
{
  set.seed(99)
  t<-yai(x=x, y=y_TH, method="randomForest", k=i)
  timp<-impute(t, method="dstWeighted", k=i, observed=T)
  timp$err<-timp$TopHeight - timp$TopHeight.o
  rmse<-sqrt(mean(timp$err^2))
  results99<-rbind(results99, c(i, rmse))
}
colnames(results99)<- c('k', 'rmse')
results99<-results99[with(results99, order(rmse)), ] #Order by RMSE
print(paste("Optimal value for k is ", results99 [1,1], sep=""))

#Model with optimal k
#TH.imp<-yai(x=x, y=y_TH, method="randomForest", k=results99[1,1], ntree=100000)
#TH.imp.out<-impute(TH.imp, vars="TopHeight", observed=TRUE, method='dstWeighted')


#####################################################
#Try an imputation using the 5 top predictors from each model
#######################################################

TH.varbs<-yaiRFsummary(TH.imp, nTop=3)
BA.varbs<-yaiRFsummary(BA.imp, nTop=3)
SPH.varbs<-yaiRFsummary(SPH.imp, nTop=3)
TRV.varbs<-yaiRFsummary(TRV.imp, nTop=3)
PS.varbs<-yaiRFsummary(P.imp, nTop=3)
UP.varbs<-yaiRFsummary(UP.imp, nTop=3)
Pulp.varbs<-yaiRFsummary(Pulp.imp, nTop=3)


TH.names<-names(as.data.frame(TH.varbs$scaledImportance))
BA.names<-names(as.data.frame(BA.varbs$scaledImportance))
SPH.names<-names(as.data.frame(SPH.varbs$scaledImportance))
TRV.names<-names(as.data.frame(TRV.varbs$scaledImportance))
PS.names<-names(as.data.frame(PS.varbs$scaledImportance))
UP.names<-names(as.data.frame(UP.varbs$scaledImportance))
Pulp.names<-names(as.data.frame(Pulp.varbs$scaledImportance))

varbs<-rbind(TH.names, BA.names)
varbs<-rbind(varbs, SPH.names)
varbs<-rbind(varbs, TRV.names)
varbs<-rbind(varbs, PS.names)
varbs<-rbind(varbs, UP.names)
varbs<-rbind(varbs, Pulp.names)

varbs<-melt(varbs)

#This will add only the top 5 from each separate imputation into x
l<-varbs$value
l<-match(l, names(final))
x=final[,l] # output from the top 5 excercise
write.csv(x, paste(output.dir, 'x.csv', sep=""))
