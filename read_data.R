### Read in Tairua data from csv files supplied by Jonathan
### d4 = N+n row data.frame with pixel information, including identifier, location (x,y), forest,cpt,stand and estyear
### d3 = n row data.frame with reference plot information, including identifier and response variable values
### d1 = N+n * n data.frame. Each row contains n reference plot identifiers for target pixels and also for reference plots  ordered by distance in covariate space from the pixel/plot identified by the row name
### d2 = N+n * n data.frame containing the distances in co-variate space corresponding to the nearest neighbours in d1
### metrics is N+N * q data.frame containing q LiDAR metrics for each of the N target pixels and n reference plots
### d1-d4, metrics use plot/pixel identifier as rowname
### reference.ids is character vector of reference plot identifiers
### target.ids is character vector of target plot identifiers
### extra.reference.ids is character vector of purposively selected reference plot identifiers
library(nlme)
dataFolder <- data.dir
metricsFolder <- data.dir
referenceMetrics.name <- "" # If reference metrics are in separate file
### Locations and compartment id of each pixel
### This assumes that each location has only one pixel
### The data includes both target and reference pixels
### It is not complete with respect to identifiers or guaranteed to be ordered
d4 <- read.csv(paste(dataFolder,"VarianceHeader.csv",sep=""))
names(d4) <- tolower(make.names(names(d4),allow_=F))
rownames(d4) <- d4$identifier
### Dave Pont only sampled stands between ages 11 and 35; estyear 1977-2000
### But there are 5 plots in 2001 stands
###d4 <- subset(d4, identifier > 400)

### It seems that the forest name differs depending on whether this is a target or a reference.  Assume all in same forest
d4$compartment <- as.factor(d4$compartment)
d4$stand <- as.factor(d4$stand)
d4$stand.id <-  getGroups(d4,~compartment/stand,level=2)
d4 <- d4[ order(d4$identifier),]

### Swap x,y co-ordinates where necessary
### Easting should be less than Northing
x <- with(d4,ifelse(center.x < center.y, center.x, center.y))
y <- with(d4,ifelse(center.x < center.y, center.y, center.x))
d4 <- within(d4, {
  center.x <- x 
  center.y <- y 
})
rm(x,y)



### Reference pixels have ids 1-99
### These are the reference pixel ids used to index the distances
### This table contains multiple pixels at each location because of the way it was generated
d1 <- read.table(paste(dataFolder,'DistanceIDS.csv',sep=""),sep=',',header=T,row.names="target",na.strings=c("NA","NULL"))
names(d1) <- tolower(names(d1))
d1 <- d1[ order(as.numeric(rownames(d1))),]


### Response variables  by reference pixel
d3 <- read.csv(paste(dataFolder,"Reference_Variance.csv",sep=""))
names(d3)[ names(d3) == "UID" ] <- "identifier"
rownames(d3) <- d3$identifier
d3 <- d3[ order(d3$identifier), ]
### Check that it is safe to use a numeric index in lieu of a character index
#stopifnot(all( as.numeric(rownames(d3)) == 1:nrow(d3)))



### Generate complete set of target pixel identifiers as character for indexing incomplete knn matrix
target.ids <- as.character( sort(intersect(as.numeric(rownames(d1)), unique(subset( d4$identifier, d4$identifier > 196 )))))
### All reference points
reference.ids <- as.character( sort( subset( d4$identifier, d4$identifier < 197 )))

### Reference points that were purposely selected; not randomly sampled.  Don't know these
extra.reference.ids <- as.character(c())

### Reference points that were part of systematic sample
gridded.reference.ids <- reference.ids[ ! reference.ids %in% extra.reference.ids ]


### Finally, subset d4 so that it only contains pixels that are in the sampling frame
d4 <- d4[ c(reference.ids,target.ids),]
d4$stand.id <- d4$stand.id[drop=T]



