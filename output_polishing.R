# **********************************************************************************************
# This Script will clip the output raster to the sampling frame and any final 
# J. Dash
# Jonathan.dash@scionresearch.com 
# March 2016
# ***********************************************************************************************

# Set libraries

library(dplyr)
library(tidyr)
library(reshape2)
library(yaImpute)
library(sp)
library(rgeos)
library(rgdal)
library(raster)
library(ggplot2)

# Read data

# Get the shape containing anciliary data
sf<-shapefile(paste(data.dir, 'SamplingFrame.shp', sep = ""))
proj4string(sf)<-proj4string

# This won't be needed eventually
alt.imp.out.targ<-stack(paste(output.dir, ' Matakana_responses.tif', sep=""))
plot(alt.imp.out.targ)

names(alt.imp.out.targ)<-c("BasalArea", "TopHeight", "Stocking", 
                           "TotalRecoverableVolume", "P_Sawlog", "UP_Sawlog", "Pulp")

## crop and mask
r2 <- crop(alt.imp.out.targ, extent(sf))
r3 <- mask(r2, sf)
plot(r3)


writeRaster(r3, filename = paste(output.dir,names(r3)), bylayer=TRUE, format="GTiff", overwrite=TRUE)

