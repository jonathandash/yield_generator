#*******************************************************************
# This module will finalise the target dataset for the matakana inventory,
# this includes pixel level tree counts etc. A finalised dataframe suitable  
# for imputation with yaImpute is produced. This requires that matakana_get_reference.R
# is executed to generate the reference dataset.
# J.Dash (jonathan.dash@scionresearch.com) Jan 2016
#*******************************************************************

library(dplyr)
library(sp)
library(rgeos)

#target.lid<-subset(target.lid, select=-c(cov_gap, dns_gap))
#target.lid<-read.csv('mata_lidar.csv',  stringsAsFactors = F, na.strings = "-"))
target.lid<-read.csv(paste(data.dir,'mata_lidar.csv', sep=""),  stringsAsFactors = F, na.strings = "-")
#Make the target lidar a raster
coordinates(target.lid)<- ~x+y
gridded(target.lid)=TRUE



###This is for the EABA stuff
#Get the extent of the the target raster
ex<-extent(target.lid)

#Create an empty raster based on the extent of the target
tr.counts.ras<-raster(ext=ex, resolution=25)

# Rasterize the Individual tree data
# Tree counts first
loc.max.tr = loc.max[c("Tree.no.")] #Limit the data for a tree count
tr.counts.count<-rasterize(loc.max.tr, tr.counts.ras, fun='count') #Count the entries in each cell 
names(tr.counts.count)<-c("ID", "p.count")
plot(tr.counts.count$p.count) #Plot the tree counts per pixel for a quick logic check
print("This figure shows an estimate of the number of trees per 25 m pixel")

# Max crown height now
loc.max.mean = loc.max[c("Height", "Max.crown.width", "Ht.to.crown.base")]
tr.counts.max<-rasterize(loc.max.mean, tr.counts.ras, fun='max')
names(tr.counts.max)<-c("ID", "max.chm.ht", "max.cr.width", "max.Ht.cr.base")
#plot(tr.counts.max)

# Min crown base now
tr.counts.min<-rasterize(loc.max.mean, tr.counts.ras, fun='min')
names(tr.counts.min)<-c("ID", "min.chm.ht", "min.cr.width", "min.Ht.cr.base")
#plot(tr.counts.min)

# Mean
loc.max.mean = loc.max[c("Height", "Max.crown.width", "Ht.to.crown.base")] #Limit the data for max crown width
tr.counts.mean<-rasterize(loc.max.mean, tr.counts.ras, fun=mean)
names(tr.counts.mean)<-c("ID", "mean.chm.ht", "mean.cr.width", "mean.Ht.cr.base")
#plot(tr.counts.mean)

# standard deviation
tr.counts.std<-rasterize(loc.max.mean, tr.counts.ras, fun=sd)
names(tr.counts.std)<-c("ID", "std.chm.ht", "std.cr.width", "std.Ht.cr.base")

# variance
tr.counts.var<-rasterize(loc.max.mean, tr.counts.ras, fun=var)
names(tr.counts.var)<-c("ID", "var.chm.ht", "var.cr.width", "var.Ht.cr.base")

#Put all the ITA results in a stack
targ.ITA<-stack(tr.counts.count$p.count, 
                tr.counts.min$min.chm.ht, tr.counts.min$min.cr.width, tr.counts.min$min.Ht.cr.base,
                tr.counts.mean$mean.chm.ht, tr.counts.mean$mean.cr.width, tr.counts.mean$mean.Ht.cr.base,
                tr.counts.max$max.chm.ht, tr.counts.max$max.cr.width, tr.counts.max$max.Ht.cr.base,
                tr.counts.std$std.chm.ht, tr.counts.std$std.cr.width, tr.counts.std$std.Ht.cr.base,
                tr.counts.var$var.chm.ht, tr.counts.var$var.cr.width, tr.counts.var$var.Ht.cr.base)


###Now get the attributes from the shapefile
for (i in 1:(nrow(shp@data)))
{
  shp@data$ID[i] = (i+1000)
}

standattribute = shp@data[c("ID", "NZ_Site", "Avg_Elev", "Avg_Slp_De", "Avg_Aspect", "Cur_silvi_", "Age", "StandKey", "Species", "Forest", "Compartmen", "Stand", "Yearplant")]
r1<-rasterize(shp, tr.counts.ras, field="ID")

rp = rasterToPoints(r1,spatial=TRUE) 
dfrp<-as.data.frame(rp)

target.lid.df<-as.data.frame(target.lid)

target.lid.df.aux<-merge(target.lid.df, dfrp,by = c("x","y"))
standattribute$newID = standattribute$ID

target_extended = merge(target.lid.df.aux, standattribute,by.x = "layer",by.y = "ID")

#Convert ITA stack to df
targ.ITA<-as.data.frame(targ.ITA, xy=T)
target_extended<-merge(target_extended, targ.ITA,by = c("x","y"))

target_extended[is.na(target_extended)]<-0

#Add in Null spaces for the response variables
target_extended$TopHeight <- NA
target_extended$BasalArea <- NA
target_extended$TotalRecoverableVolume <- NA
target_extended$Stocking <- NA
target_extended$P_Sawlog <- NA
target_extended$UP_Sawlog <- NA
target_extended$Pulp<- NA


target_extended$p.stocking<-target_extended$p.count/0.06
#Calculate interactions for target

target_extended[target_extended==0]<-0.0000001

#Age interaction
target_extended$ag.avg <- target_extended$avg*target_extended$Age
target_extended$ag.kur <- target_extended$kur*target_extended$Age
target_extended$ag.p30 <- target_extended$p30*target_extended$Age
target_extended$ag.p95 <- target_extended$p95*target_extended$Age
target_extended$ag.b20 <- target_extended$b20*target_extended$Age
target_extended$ag.b90 <- target_extended$b90*target_extended$Age
target_extended$ag.c01 <- target_extended$c01*target_extended$Age
target_extended$ag.c05 <- target_extended$c05*target_extended$Age
target_extended$ag.dns <- target_extended$dns*target_extended$Age
target_extended$ag.Avg_Aspect <- target_extended$Avg_Aspect*target_extended$Age
target_extended$ag.var.cr.width <- target_extended$var.cr.width*target_extended$Age
target_extended$ag.mean.Ht.cr.base <- target_extended$mean.Ht.cr.base*target_extended$Age
target_extended$ag.var.Ht.cr.base <- target_extended$var.Ht.cr.base*target_extended$Age
target_extended$ag.std.chm.ht <- target_extended$std.chm.ht*target_extended$Age
target_extended$ag.Age <- target_extended$Age*target_extended$Age
target_extended$ag.qav <- target_extended$qav*target_extended$Age
target_extended$ag.p05 <- target_extended$p05*target_extended$Age
target_extended$ag.p50 <- target_extended$p50*target_extended$Age
target_extended$ag.p99 <- target_extended$p99*target_extended$Age
target_extended$ag.b30 <- target_extended$b30*target_extended$Age
target_extended$ag.b95 <- target_extended$b95*target_extended$Age
target_extended$ag.c02 <- target_extended$c02*target_extended$Age
target_extended$ag.c06 <- target_extended$c06*target_extended$Age
target_extended$ag.NZ_Site <- target_extended$NZ_Site*target_extended$Age
target_extended$ag.std.cr.width <- target_extended$std.cr.width*target_extended$Age
target_extended$ag.min.Ht.cr.base <- target_extended$min.Ht.cr.base*target_extended$Age
target_extended$ag.mean.chm.ht <- target_extended$mean.chm.ht*target_extended$Age
target_extended$ag.var.chm.ht <- target_extended$var.chm.ht*target_extended$Age
target_extended$ag.min <- target_extended$min*target_extended$Age
target_extended$ag.std <- target_extended$std*target_extended$Age
target_extended$ag.p10 <- target_extended$p10*target_extended$Age
target_extended$ag.p70 <- target_extended$p70*target_extended$Age
target_extended$ag.b05 <- target_extended$b05*target_extended$Age
target_extended$ag.b50 <- target_extended$b50*target_extended$Age
target_extended$ag.b99 <- target_extended$b99*target_extended$Age
target_extended$ag.c03 <- target_extended$c03*target_extended$Age
target_extended$ag.c07 <- target_extended$c07*target_extended$Age
target_extended$ag.Avg_Elev <- target_extended$Avg_Elev*target_extended$Age
target_extended$ag.p.count <- target_extended$p.count*target_extended$Age
target_extended$ag.min.cr.width <- target_extended$min.cr.width*target_extended$Age
target_extended$ag.max.Ht.cr.base <- target_extended$max.Ht.cr.base*target_extended$Age
target_extended$ag.min.chm.ht <- target_extended$min.chm.ht*target_extended$Age
target_extended$ag.p.stocking <- target_extended$p.stocking*target_extended$Age
target_extended$ag.max <- target_extended$max*target_extended$Age
target_extended$ag.ske <- target_extended$ske*target_extended$Age
target_extended$ag.p20 <- target_extended$p20*target_extended$Age
target_extended$ag.p90 <- target_extended$p90*target_extended$Age
target_extended$ag.b10 <- target_extended$b10*target_extended$Age
target_extended$ag.b70 <- target_extended$b70*target_extended$Age
target_extended$ag.c00 <- target_extended$c00*target_extended$Age
target_extended$ag.c04 <- target_extended$c04*target_extended$Age
target_extended$ag.cov <- target_extended$cov*target_extended$Age
target_extended$ag.Avg_Slp_De <- target_extended$Avg_Slp_De*target_extended$Age
target_extended$ag.mean.cr.width <- target_extended$mean.cr.width*target_extended$Age
target_extended$ag.max.cr.width <- target_extended$max.cr.width*target_extended$Age
target_extended$ag.std.Ht.cr.base <- target_extended$std.Ht.cr.base*target_extended$Age
target_extended$ag.max.chm.ht <- target_extended$max.chm.ht*target_extended$Age


#Squares
target_extended$avg.2 <- target_extended$avg^2
target_extended$kur.2 <- target_extended$kur^2
target_extended$p30.2 <- target_extended$p30^2
target_extended$p95.2 <- target_extended$p95^2
target_extended$b20.2 <- target_extended$b20^2
target_extended$b90.2 <- target_extended$b90^2
target_extended$c01.2 <- target_extended$c01^2
target_extended$c05.2 <- target_extended$c05^2
target_extended$dns.2 <- target_extended$dns^2
target_extended$Avg_Aspect.2 <- target_extended$Avg_Aspect^2
target_extended$var.cr.width.2 <- target_extended$var.cr.width^2
target_extended$mean.Ht.cr.base.2 <- target_extended$mean.Ht.cr.base^2
target_extended$var.Ht.cr.base.2 <- target_extended$var.Ht.cr.base^2
target_extended$std.chm.ht.2 <- target_extended$std.chm.ht^2
target_extended$Age.2 <- target_extended$Age^2
target_extended$qav.2 <- target_extended$qav^2
target_extended$p05.2 <- target_extended$p05^2
target_extended$p50.2 <- target_extended$p50^2
target_extended$p99.2 <- target_extended$p99^2
target_extended$b30.2 <- target_extended$b30^2
target_extended$b95.2 <- target_extended$b95^2
target_extended$c02.2 <- target_extended$c02^2
target_extended$c06.2 <- target_extended$c06^2
target_extended$NZ_Site.2 <- target_extended$NZ_Site^2
target_extended$std.cr.width.2 <- target_extended$std.cr.width^2
target_extended$min.Ht.cr.base.2 <- target_extended$min.Ht.cr.base^2
target_extended$mean.chm.ht.2 <- target_extended$mean.chm.ht^2
target_extended$var.chm.ht.2 <- target_extended$var.chm.ht^2
target_extended$min.2 <- target_extended$min^2
target_extended$std.2 <- target_extended$std^2
target_extended$p10.2 <- target_extended$p10^2
target_extended$p70.2 <- target_extended$p70^2
target_extended$b05.2 <- target_extended$b05^2
target_extended$b50.2 <- target_extended$b50^2
target_extended$b99.2 <- target_extended$b99^2
target_extended$c03.2 <- target_extended$c03^2
target_extended$c07.2 <- target_extended$c07^2
target_extended$Avg_Elev.2 <- target_extended$Avg_Elev^2
target_extended$p.count.2 <- target_extended$p.count^2
target_extended$min.cr.width.2 <- target_extended$min.cr.width^2
target_extended$max.Ht.cr.base.2 <- target_extended$max.Ht.cr.base^2
target_extended$min.chm.ht.2 <- target_extended$min.chm.ht^2
target_extended$p.stocking.2 <- target_extended$p.stocking^2
target_extended$max.2 <- target_extended$max^2
target_extended$ske.2 <- target_extended$ske^2
target_extended$p20.2 <- target_extended$p20^2
target_extended$p90.2 <- target_extended$p90^2
target_extended$b10.2 <- target_extended$b10^2
target_extended$b70.2 <- target_extended$b70^2
target_extended$c00.2 <- target_extended$c00^2
target_extended$c04.2 <- target_extended$c04^2
target_extended$cov.2 <- target_extended$cov^2
target_extended$Avg_Slp_De.2 <- target_extended$Avg_Slp_De^2
target_extended$mean.cr.width.2 <- target_extended$mean.cr.width^2
target_extended$max.cr.width.2 <- target_extended$max.cr.width^2
target_extended$std.Ht.cr.base.2 <- target_extended$std.Ht.cr.base^2
target_extended$max.chm.ht.2 <- target_extended$max.chm.ht^2


#sqrt
target_extended$avg.05 <- sqrt(abs(target_extended$avg))
target_extended$kur.05 <- sqrt(abs(target_extended$kur))
target_extended$p30.05 <- sqrt(abs(target_extended$p30))
target_extended$p95.05 <- sqrt(abs(target_extended$p95))
target_extended$b20.05 <- sqrt(abs(target_extended$b20))
target_extended$b90.05 <- sqrt(abs(target_extended$b90))
target_extended$c01.05 <- sqrt(abs(target_extended$c01))
target_extended$c05.05 <- sqrt(abs(target_extended$c05))
target_extended$dns.05 <- sqrt(abs(target_extended$dns))
target_extended$Avg_Aspect.05 <- sqrt(abs(target_extended$Avg_Aspect))
target_extended$var.cr.width.05 <- sqrt(abs(target_extended$var.cr.width))
target_extended$mean.Ht.cr.base.05 <- sqrt(abs(target_extended$mean.Ht.cr.base))
target_extended$var.Ht.cr.base.05 <- sqrt(abs(target_extended$var.Ht.cr.base))
target_extended$std.chm.ht.05 <- sqrt(abs(target_extended$std.chm.ht))
target_extended$Age.05 <- sqrt(abs(target_extended$Age))
target_extended$qav.05 <- sqrt(abs(target_extended$qav))
target_extended$p05.05 <- sqrt(abs(target_extended$p05))
target_extended$p50.05 <- sqrt(abs(target_extended$p50))
target_extended$p99.05 <- sqrt(abs(target_extended$p99))
target_extended$b30.05 <- sqrt(abs(target_extended$b30))
target_extended$b95.05 <- sqrt(abs(target_extended$b95))
target_extended$c02.05 <- sqrt(abs(target_extended$c02))
target_extended$c06.05 <- sqrt(abs(target_extended$c06))
target_extended$NZ_Site.05 <- sqrt(abs(target_extended$NZ_Site))
target_extended$std.cr.width.05 <- sqrt(abs(target_extended$std.cr.width))
target_extended$min.Ht.cr.base.05 <- sqrt(abs(target_extended$min.Ht.cr.base))
target_extended$mean.chm.ht.05 <- sqrt(abs(target_extended$mean.chm.ht))
target_extended$var.chm.ht.05 <- sqrt(abs(target_extended$var.chm.ht))
target_extended$min.05 <- sqrt(abs(target_extended$min))
target_extended$std.05 <- sqrt(abs(target_extended$std))
target_extended$p10.05 <- sqrt(abs(target_extended$p10))
target_extended$p70.05 <- sqrt(abs(target_extended$p70))
target_extended$b05.05 <- sqrt(abs(target_extended$b05))
target_extended$b50.05 <- sqrt(abs(target_extended$b50))
target_extended$b99.05 <- sqrt(abs(target_extended$b99))
target_extended$c03.05 <- sqrt(abs(target_extended$c03))
target_extended$c07.05 <- sqrt(abs(target_extended$c07))
target_extended$Avg_Elev.05 <- sqrt(abs(target_extended$Avg_Elev))
target_extended$p.count.05 <- sqrt(abs(target_extended$p.count))
target_extended$min.cr.width.05 <- sqrt(abs(target_extended$min.cr.width))
target_extended$max.Ht.cr.base.05 <- sqrt(abs(target_extended$max.Ht.cr.base))
target_extended$min.chm.ht.05 <- sqrt(abs(target_extended$min.chm.ht))
target_extended$p.stocking.05 <- sqrt(abs(target_extended$p.stocking))
target_extended$max.05 <- sqrt(abs(target_extended$max))
target_extended$ske.05 <- sqrt(abs(target_extended$ske))
target_extended$p20.05 <- sqrt(abs(target_extended$p20))
target_extended$p90.05 <- sqrt(abs(target_extended$p90))
target_extended$b10.05 <- sqrt(abs(target_extended$b10))
target_extended$b70.05 <- sqrt(abs(target_extended$b70))
target_extended$c00.05 <- sqrt(abs(target_extended$c00))
target_extended$c04.05 <- sqrt(abs(target_extended$c04))
target_extended$cov.05 <- sqrt(abs(target_extended$cov))
target_extended$Avg_Slp_De.05 <- sqrt(abs(target_extended$Avg_Slp_De))
target_extended$mean.cr.width.05 <- sqrt(abs(target_extended$mean.cr.width))
target_extended$max.cr.width.05 <- sqrt(abs(target_extended$max.cr.width))
target_extended$std.Ht.cr.base.05 <- sqrt(abs(target_extended$std.Ht.cr.base))
target_extended$max.chm.ht.05 <- sqrt(abs(target_extended$max.chm.ht))


#log
target_extended$logavg <- log(abs(target_extended$avg))
target_extended$logkur <- log(abs(target_extended$kur))
target_extended$logp30 <- log(abs(target_extended$p30))
target_extended$logp95 <- log(abs(target_extended$p95))
target_extended$logb20 <- log(abs(target_extended$b20))
target_extended$logb90 <- log(abs(target_extended$b90))
target_extended$logc01 <- log(abs(target_extended$c01))
target_extended$logc05 <- log(abs(target_extended$c05))
target_extended$logdns <- log(abs(target_extended$dns))
target_extended$logAvg_Aspect <- log(abs(target_extended$Avg_Aspect))
target_extended$logvar.cr.width <- log(abs(target_extended$var.cr.width))
target_extended$logmean.Ht.cr.base <- log(abs(target_extended$mean.Ht.cr.base))
target_extended$logvar.Ht.cr.base <- log(abs(target_extended$var.Ht.cr.base))
target_extended$logstd.chm.ht <- log(abs(target_extended$std.chm.ht))
target_extended$logAge <- log(abs(target_extended$Age))
target_extended$logqav <- log(abs(target_extended$qav))
target_extended$logp05 <- log(abs(target_extended$p05))
target_extended$logp50 <- log(abs(target_extended$p50))
target_extended$logp99 <- log(abs(target_extended$p99))
target_extended$logb30 <- log(abs(target_extended$b30))
target_extended$logb95 <- log(abs(target_extended$b95))
target_extended$logc02 <- log(abs(target_extended$c02))
target_extended$logc06 <- log(abs(target_extended$c06))
target_extended$logNZ_Site <- log(abs(target_extended$NZ_Site))
target_extended$logstd.cr.width <- log(abs(target_extended$std.cr.width))
target_extended$logmin.Ht.cr.base <- log(abs(target_extended$min.Ht.cr.base))
target_extended$logmean.chm.ht <- log(abs(target_extended$mean.chm.ht))
target_extended$logvar.chm.ht <- log(abs(target_extended$var.chm.ht))
target_extended$logmin <- log(abs(target_extended$min))
target_extended$logstd <- log(abs(target_extended$std))
target_extended$logp10 <- log(abs(target_extended$p10))
target_extended$logp70 <- log(abs(target_extended$p70))
target_extended$logb05 <- log(abs(target_extended$b05))
target_extended$logb50 <- log(abs(target_extended$b50))
target_extended$logb99 <- log(abs(target_extended$b99))
target_extended$logc03 <- log(abs(target_extended$c03))
target_extended$logc07 <- log(abs(target_extended$c07))
target_extended$logAvg_Elev <- log(abs(target_extended$Avg_Elev))
target_extended$logp.count <- log(abs(target_extended$p.count))
target_extended$logmin.cr.width <- log(abs(target_extended$min.cr.width))
target_extended$logmax.Ht.cr.base <- log(abs(target_extended$max.Ht.cr.base))
target_extended$logmin.chm.ht <- log(abs(target_extended$min.chm.ht))
target_extended$logp.stocking <- log(abs(target_extended$p.stocking))
target_extended$logmax <- log(abs(target_extended$max))
target_extended$logske <- log(abs(target_extended$ske))
target_extended$logp20 <- log(abs(target_extended$p20))
target_extended$logp90 <- log(abs(target_extended$p90))
target_extended$logb10 <- log(abs(target_extended$b10))
target_extended$logb70 <- log(abs(target_extended$b70))
target_extended$logc00 <- log(abs(target_extended$c00))
target_extended$logc04 <- log(abs(target_extended$c04))
target_extended$logcov <- log(abs(target_extended$cov))
target_extended$logAvg_Slp_De <- log(abs(target_extended$Avg_Slp_De))
target_extended$logmean.cr.width <- log(abs(target_extended$mean.cr.width))
target_extended$logmax.cr.width <- log(abs(target_extended$max.cr.width))
target_extended$logstd.Ht.cr.base <- log(abs(target_extended$std.Ht.cr.base))
target_extended$logmax.chm.ht <- log(abs(target_extended$max.chm.ht))



#Now get rid of excess columns
target.index<-target_extended[,c('x', 'y', 'StandKey', 'Species', "Forest", "Compartmen", "Stand", "Yearplant")]
target_extended<-subset(target_extended, select=-c(layer, index, newID, StandKey, Species, Forest, Compartmen, Stand, Yearplant)) #ditch these for now


#Add In Source column
ref$Source<-"Ref"
target_extended$Source<-"Target"
paste('Column number difference is', ncol(ref) - ncol(target_extended), '... If column difference is 0 proceed')


# Target is now possibly ready for an imputation
#str(target_extended)

final<-as.data.frame(rbind(ref, target_extended))
final$Source<-as.factor(final$Source)
final$UID<-row.names(final)

#str(final)
