#******************************************************************************************
# Generate yield tables
# Jonathan Dash (jonathan.dash@scionresearch.com)
# Jan 2016
#******************************************************************************************


#Get the distances
dists<-as.data.frame(alt.imp$neiIdsTrgs)

#RENAME these ugly column names

#Add in x and y
NN_rast<-cbind(target.index, dists)

NN_melt<-melt(NN_rast, id.vars=c("x", "y", "StandKey", "Species"))

#Make the output a spatial pixel dataframe
coordinates(NN_rast)<-~x+y
gridded(NN_rast)<-TRUE
#NN_rast<-stack(NN_rast)
#plot(NN_rast)


####################################################################
#Now use the NN raster to build yield tables
#******************************************************************

#Read in future yields for reference plots
ref.stand.summary.yields<-read.csv(paste(data.dir, 'stand_summary_yields.csv', sep="")) #This is the stand summary table from ytgen

#head(ref.stand.summary)
ref.stand.summary.yields<-subset(ref.stand.summary.yields, select=c(PlotName, PeriodNumber, Age, BasalArea, TopHeight, Stocking, TotalRecoverableVolume))
names(ref.stand.summary.yields)[names(ref.stand.summary.yields) == 'PlotName'] <- 'Plot_ID' #Rename the Plot_ID column



ref.log.grade.yields<-read.csv(paste(data.dir, 'loggrade_pivot_lpm1_yields.csv', sep="")) #This table contains a one line summary for each plot detailing the 
#log product volumes required for prediction
names(ref.log.grade.yields)[names(ref.log.grade.yields) == 'PlotName'] <- 'Plot_ID' #Rename the Plot_ID column
ref.log.grade.yields<-subset(ref.log.grade.yields, select=-c(PopulationName))

#Create ref yields by merging stand and grade pivot
ref.yields<-merge(ref.stand.summary.yields, ref.log.grade.yields, by=c("Plot_ID", "PeriodNumber"))

#create an index from the imputation ref plots
ref.plots$UID<-row.names(ref.plots)
ref.index<-data.frame(UID = ref.plots$UID, Plot_ID = ref.plots$Plot_ID)

ref.yt<-merge(NN_melt, ref.index, by.x="value", by.y="UID")
ref.yt<-merge(ref.yt, ref.yields, by="Plot_ID")

#Make Yield Tables
ref.yt.out<-subset(ref.yt, ref.yt$variable %in% c("Id.k1", "id.k3", "id.k2"))
ref.yt.out.summed<- ref.yt.out %>% group_by(StandKey, Species, PeriodNumber) %>% summarise(Age = mean(Age), BasalArea=mean(BasalArea),
                                                                                           TopHeight=mean(TopHeight), TotalRecoverableVolume=mean(TotalRecoverableVolume),
                                                                                           Stocking=mean(Stocking), yP1PR=mean(yP1PR),
                                                                                           YSaw1PR=mean(ySaw1PR), ySaw2PR=mean(ySaw2PR),
                                                                                           ySaw3PR=mean(ySaw3PR), yPulp=mean(yPulp))


ref.yt.out.summed<- subset(ref.yt.out.summed, ref.yt.out.summed$Species == 'P.RAD')
ref.yt.out.summed$Year<-ref.yt.out.summed$PeriodNumber + 2014

write.csv(ref.yt.out.summed, paste(output.dir, 'yield_tables_v1.csv', sep=""), row.names = FALSE)







