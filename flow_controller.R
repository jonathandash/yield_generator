#***************************************************************************
# Flow Controller - Matakana Yield Tables
# This is the controller for the Matakana yield tables
# Jonathan Dash January 2016 (jonathan.dash@scionresearch.com)
#***************************************************************************

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


#String relating to nztm
proj4string <- "+proj=tmerc +lat_0=0 +lon_0=173 +k=0.9996 +x_0=1600000 +y_0=10000000 +ellps=GRS80 +units=m +no_defs"


# Set up work directories
# This should be automated to the location of the controller
setwd("D:\\MatakanaIsland\\R scripts\\Finalised")

# set data location - Assumes all input data are in a directory called "data"
data.dir<-"./data/"

# Create an output directory called outputs to store all outputs
dir.create("./outputs/", showWarnings = TRUE)
output.dir<- "./outputs/"


#*******************************************************************************
# Set variables for imputation
#*******************************************************************************

# Select a value of k (This can be automated)
k <- 3

# Select a method for defining nearest neigbours
method<-"randomForest"

# Select imputation methods
imp<-"dstWeighted"

# Select responses
responses<-c("BasalArea", "TopHeight", "Stocking", "TotalRecoverableVolume", "P_Sawlog", "UP_Sawlog", "Pulp")

# Variable selection method

# Number of variables

# Call other R modules
source('matakana_get_reference.R') #Get reference data
source('matakana_get_target.R') #Get target data
source('matakana_varb_sel_top5.R') #Variable selection
source('matakana_knn_imputation.R') #Modelling
source('matakana_yield_tables.R') #Make yield tables
source('matakana_generate_se_inputs.R') #Get inputs for se calculation
source('matakana_SE.R') # Calculate stand level standard error 
source('matakana_output_polishing.R')




