### Read old version of estyear by stand.id
### This is because it isn't included in updates of data
get.estyear <- function()
  {
    dataFolder <- "./data/"
    d <- read.csv(paste(dataFolder,"AncientVarianceHeader.csv",sep=""))
    names(d) <- tolower(names(d))
    d$compartment <- as.factor(d$compartment)
    d$stand.id <-  getGroups(d,~forest/compartment/stand,level=3)

    with(d,tapply(estyear,stand.id,max))

  }
