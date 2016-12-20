#******************************************************************************************
# Carry out an imputation and output results
# Jonathan Dash (jonathan.dash@scionresearch.com)
# Jan 2016
#******************************************************************************************

# These should be specified in the flow controller
#responses<-c("BasalArea", "TopHeight", "Stocking", "TotalRecoverableVolume", "P_Sawlog", "UP_Sawlog", "Pulp")

y<-as.data.frame(subset(final, final$Source=="Ref"))


# These should be specified in the flow controller
y<-y[c("BasalArea", "TopHeight", "Stocking", "TotalRecoverableVolume", "P_Sawlog", "UP_Sawlog", "Pulp")]

#Carry out an imputation using the above
alt.imp<-yai(x=x, y=y, method="randomForest", k=k)
plot(alt.imp, vars = responses)
z<-rmsd(alt.imp, vars = responses, scale=T)
#grmsd(alt.imp, vars = responses)

tiff('D:\\MatakanaIsland\\final_analysis_yieldtables\\SamplingError\\Paper\\rmsd_s.tif', w=12, h=10, units="cm", res=300)
#par(mar = c(0, 1.5, 0, 0), oma = c(4, 4, 3, 1), tcl = 0.35, mgp = c(2, 0.4, 0))
barplot(z$rmsdS, names = c("BA", "TH", "Sph", "TRV", "PS", "UPS", "Pulp"), space = 0.5, las=1, ylab="Scaled RMSD", ylim=c(0,1.2))
box()
dev.off()

#Impute the responses
alt.imp.out<-impute(alt.imp, vars=responses, observed=TRUE, method='dstWeighted')


#Filter out the reference plots
alt.imp.out.targ<-(alt.imp.out[!complete.cases(alt.imp.out), ])

#Filter out the target plots
alt.imp.out.ref<-(alt.imp.out[complete.cases(alt.imp.out), ])

# Plot the imputation for inclusion in the paper
tiff('D:\\MatakanaIsland\\final_analysis_yieldtables\\SamplingError\\Paper\\impute_plots.tif', w=21 , h=15, units="cm", res=300)
layout(matrix(1:6, nrow = 2))
par(mar = c(3, 3.5, 0, 0), oma = c(4, 4, 3, 1), tcl = 0.35, mgp = c(2, 0.4, 0))
plot(TotalRecoverableVolume~ TotalRecoverableVolume.o, data = alt.imp.out.ref, las=1,
     ylab=expression(Imputed~total~recoverable~volume~m^3~ha^-1),
     xlab=expression(Observed~total~recoverable~volume~m^3~ha^-1))
abline(0, 1)
plot(TopHeight~ TopHeight.o, data = alt.imp.out.ref, las=1,
     ylab=expression(Imputed~top~height~m),
     xlab=expression(Observed~top~height~m))
abline(0,1)
plot(BasalArea~ BasalArea.o, data = alt.imp.out.ref, las=1,
     ylab=expression(Imputed~basal~area~m^2~ha^-1),
     xlab=expression(Observed~basal~area~m^2~ha^-1))
abline(0,1)
plot(Stocking~ Stocking.o, data = alt.imp.out.ref, las=1,
     ylab=expression(Imputed~stems~per~hectare),
     xlab=expression(Observed~stems~per~hectare))
abline(0,1)
plot(UP_Sawlog~ UP_Sawlog.o, data = alt.imp.out.ref, las=1,
     ylab=expression(Imputed~unpruned~sawlog~volume~m^3~ha^-1),
     xlab=expression(Observed~unpruned~sawlog~volume~m^3~ha^-1))
abline(0,1)
plot(P_Sawlog~ P_Sawlog.o, data = alt.imp.out.ref, las=1,
     ylab=expression(Imputed~pruned~sawlog~volume~m^3~ha^-1),
     xlab=expression(Observed~pruned~sawlog~volume~m^3~ha^-1))
abline(0,1)
dev.off()

#Get cooordinates from the target
targ.coords<-target_extended[c('x', 'y')]

#Add the corrdinates in to the imputation output
alt.imp.out.targ<-cbind(targ.coords, alt.imp.out.targ)

#Get rid of the unwanted fields
alt.imp.out.targ<-alt.imp.out.targ[c('x', 'y', "BasalArea", "TopHeight", "Stocking", 
                                     "TotalRecoverableVolume", "P_Sawlog", "UP_Sawlog", "Pulp")]

#Make the output a raster stack
coordinates(alt.imp.out.targ)<-~x+y
gridded(alt.imp.out.targ)<-TRUE
alt.imp.out.targ<-stack(alt.imp.out.targ)
plot(alt.imp.out.targ)

#Save the stack so can use it in the paper analysis
writeRaster(alt.imp.out.targ, paste(output.dir, "Matakana_responses.tif"))

# alt.surface.sums<-extract(alt.imp.out.targ, shp, method='simple', fun=mean, df=T, na.rm=T) #density
# stand.output<-cbind(alt.surface.sums, shp)
# stand.output<-stand.output[c("Compartmen", "Stand" , "StandKey", "Species", "NetArea", "Stdkey", "Cur_silvi_", 
#                              "Age", "BasalArea", "TopHeight", "Stocking", "TotalRecoverableVolume", "P_Sawlog", "UP_Sawlog", "Pulp")]
# stand.output<-subset(stand.output, stand.output$Species=="P.RAD")
# stand.output<-subset(stand.output, stand.output$Age>9)
# 
# write.csv(stand.output, "Yields_MeasurementDate_ALS_Only.csv", row.names=FALSE)

