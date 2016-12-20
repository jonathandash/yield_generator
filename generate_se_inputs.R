#***********************************************************
#Produce the outputs for variance calculation 
#***********************************************************

#Carry out an imputation using the above
rfall10<-yai(x=x, y=y, method="randomForest", k=nrow(ref.plots))

####Build distance tables for variance calculation
TargetDistall<- data.frame(rfall10$neiDstTrgs)
TargetDistall$target <- rownames(TargetDistall)
RefsDistall<- data.frame(rfall10$neiDstRefs)
RefsDistall$target <- rownames(RefsDistall)
Distances<-merge(TargetDistall, RefsDistall, all=TRUE) 

TargetDistIDSall<- data.frame(rfall10$neiIdsTrg)
TargetDistIDSall$target <- rownames(TargetDistIDSall)
RefsDistIDSall<- data.frame(rfall10$neiIdsRefs)
RefsDistIDSall$target <- rownames(RefsDistIDSall)
DistanceIDS<-merge(TargetDistIDSall, RefsDistIDSall, all=TRUE)

#sqlSave(channel, Distances , tablename = "DistancesNew", append = TRUE, rownames = FALSE, colnames = FALSE, nastring = NULL)
#sqlSave(channel, DistanceIDS , tablename = "DistanceIDSNew", append = TRUE, rownames = FALSE, colnames = FALSE, nastring = NULL)

write.csv(Distances, paste(data.dir, 'Distances.csv', sep=""), row.names=FALSE)
write.csv(DistanceIDS, paste(data.dir, "DistanceIDS.csv", sep=""), row.names=FALSE)
#write.csv(f, "C:\\waiapapa_panpac\\foruse.csv")
write.csv(final, paste(data.dir, "metrics.csv", sep=""), row.names=FALSE)

Reference_Variance<-subset(final, final$Source =='Ref')
Reference_Variance<-Reference_Variance[c("UID", "BasalArea", "TopHeight", "Stocking", "TotalRecoverableVolume", "P_Sawlog", "UP_Sawlog", "Pulp", "Age")]
write.csv(Reference_Variance, paste(data.dir, "Reference_Variance.csv", sep=""), row.names=FALSE)

#Produce Variance Header
tar.head<-merge(target.index, final, by=c("x", "y"))
tar.head<-tar.head[c("x", "y", "Forest", "Compartmen", "Stand", "UID", "Yearplant")]
#Make names match the input in the example using BR script
names(tar.head)<-c("center_X", "center_Y", "Forest", "Compartment", "Stand", "Identifier", "ESTYEAR")


ref.head<-over(ref.loc, shp)
ref.head<-cbind(ref, ref.head)
ref.head$UID<-row.names(ref.head)
ref.head<-ref.head[c("x", "y", "Forest", "Compartmen", "Stand", "UID", "Yearplant")]
names(ref.head)<-c("center_X", "center_Y", "Forest", "Compartment", "Stand", "Identifier", "ESTYEAR")
var.head<-rbind(ref.head, tar.head)

write.csv(var.head, paste(data.dir, "VarianceHeader.csv", sep=""), row.names=FALSE)
