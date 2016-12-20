### This will calculate standard error for stands 
### For Matakana
library(nlme)
library(fields) #rdist
library(ggplot2)
library(tidyr)
library(stringr)
options(cores=4)
source('read_matakana_data.R')
source('knnfn.R')
source('calculate_correlation_matrix.R')

### Assumes that references ids are 1:n
stopifnot( all(reference.ids == as.character(1:length(reference.ids))))

### Reference points with Location and Y variables
rp <- merge(d4[reference.ids,],d3, by='identifier')
rp <- rp[order(rp$identifier),]
### Want to index into rp using a numeric identifier as a row index, so make sure identifiers and rows match
stopifnot(all( rp$identifier == row(rp)[,1]))
rownames(rp) <- rp$identifier

### Construct spatial correlation matrix with no spatial correlation
nr <- length(reference.ids)
rho.ij <- matrix(0,ncol=nr,nrow=nr)
diag(rho.ij) <- 1
rownames(rho.ij) <- rownames(rp)

### Stand group identifiers
g <- d4[target.ids, "stand.id"] # stand identifiers by target pixel
gn <- tapply(g,g,length)
gl <- names(gn)[!is.na(gn) & gn > 5 ] # distinct stand identifiers (levels of g) with reference pixels
#gl <- head(gl)  # Use a subset for testing --------------------Commented out by JD for full run

### N*n distance matrix for the extra correlation components in equation 15b
distance <- rdist( d4[ target.ids, c("center.x","center.y") ], d4[ reference.ids, c("center.x","center.y")] )

## Combinations of k and response variable
runs <- expand.grid(k=c(2,3,10),response.variable=c("TotalRecoverableVolume"))

#runs <- runs[1,]  ### Use a subset for testing --------------------Commented out by JD for full run

rho.ij.default <- rho.ij
### Calculate mean and variance by stand assuming spatial correlation
### Using sampling and the M1 estimator from McRoberts et al
xy <- rp[reference.ids,c("center.x","center.y")]
rownames(xy) <- reference.ids
ts <- proc.time()
results.2 <- do.call("rbind",lapply(1:nrow(runs), function(ri) {   #JD mclapply changed to lapply
  k <- as.numeric(runs[ri,"k"])
  response.var <- as.character(runs[ri,"response.variable"])
  cat("k=",k," response=",response.var,"\n")
  t1 <- as.matrix(d1[ target.ids, 1:k ])
  rv <- as.vector(rp[,response.var]) # response variable
  idk <- d1[reference.ids,1:k]
  rownames(idk) <- reference.ids
  gamma.fn <- NULL
  dil <- NULL
  l1 <- try(calculate.correlation.matrix(rv, idk,xy,1500),silent=T)
  # If we successfully modelled spatial correlation then use the generated
  # correlation matrix. Other wise use the one with no spatial correlation
    if (class(l1) == 'list' && !is.null(l1$vf) && "variogramModel" %in% class(l1$vf) && class(attr(l1$vf,'singular')) == 'logical' && attr(l1$vf,'singular') == FALSE)
    {
      rho.ij <- l1$rho.ij
      range <- l1$vf[2,"range"]
      dil <- distance
      gamma.fn <- exp2f(l1$vf) # Eq 7a semi-variogram function
      gamma.infinity <- gamma.fn(max(dil) * 2)   # gamma total from 7b
      rho.fn <- function(d) { 1 - gamma.fn(d) / gamma.infinity }
    }   else  {
      rho.ij <- rho.ij.default
      range <- 0
    }
  stopifnot(all(diag(rho.ij) == 1))

  tmp <- as.data.frame(t(simplify2array(lapply( gl, function(gi) {
    s <- g == gi
    cat("Stand=",gi,"\n")
    result <- numeric(6)
    tmp.result <- try(knn( t1[s,,drop=F], rv, rho.ij,method='sample' ,block.size=1000,sample.size=500*499/2,use.mc=F),silent=T)
    if ( class(tmp.result) == "numeric" ) { result <- tmp.result }
    result
  }))))

  within(tmp, {
    stand.id <- gl
    response <- response.var
    se <- sqrt(var)
    ple <- 2 * se / mu * 100
    range <- range
  })
}))
proc.time() - ts
tapply(results.2$ple, getGroups(results.2,~k/response,level=2),median)
finalresults=as.data.frame(results.2)
#finalresults$rep = d5$rep
#finalresults$percentage = d5$percentage
write.table(as.data.frame(finalresults),paste(output.dir,'stand_level_results.dat',sep=""),append=TRUE)


spl<-as.data.frame(str_split_fixed(finalresults$stand.id, "/", 2)) #separate cpt and stand
names(spl)<-c("Compartmen", "Stand")
finalresults<-cbind(finalresults, spl)
finalresults$Compartmen<-str_pad(finalresults$Compartmen, 3, pad = "0")
finalresults$Stand<-str_pad(finalresults$Stand, 3, pad = "0")

# ### Calculate sampling error over entire study area for several values of k
# ts <- proc.time()
# results <- lapply(c(2), function(k) {         #mclapply changed lapply
#   t1 <- as.matrix(d1[ target.ids, 1:k ])
#   yi <- rp[reference.ids,"TotalRecoverableVolume",drop=T]
#   idk <- d1[reference.ids,1:k]
#   rownames(idk) <- reference.ids
#   l1 <- try(calculate.correlation.matrix(yi, idk,xy,1500),silent=T)
#   # If we successfully modelled spatial correlation then use the generated
#   # correlation matrix. Otherwise use the one with no spatial correlation
#   if (class(l1) == 'list' && !is.null(l1$vf) && "variogramModel" %in% class(l1$vf) && class(attr(l1$vf,'singular')) == 'logical' && attr(l1$vf,'singular') == FALSE)
#     {
#       rho.ij <- l1$rho.ij
#     }   else  {
#       rho.ij <- rho.ij.default
#     }
# 
#   knn( t1, yi, rho.ij,method='sample',block.size=999999,sample.size=500*499/2,use.mc=F )
#                })
# proc.time() - ts
# #fn <- paste('total_study_area_results_2_',format(Sys.time(), "%Y-%m-%d_%H_%M"),'.dat',sep="")
# #write.table(as.data.frame(t(simplify2array(results))),fn)
# 
# finalresults1=as.data.frame(t(simplify2array(results)))
# #finalresults1$rep = d5$rep
# #finalresults1$percentage = d5$percentage
# write.table(as.data.frame(finalresults1),paste('Total_level_results.dat',sep=""),append=TRUE)
# 
# finalresults$Area<-0.0625*finalresults$N
# fin.graph<-subset(finalresults, finalresults$k==2)
# fin.graph<-subset(fin.graph, fin.graph$response=="TotalRecoverableVolume")
# #fin.graph<-subset(fin.graph, fin.graph$Area>5)
# 
# mean(fin.graph$ple)
# TotalRecoverableVolume
# 
# 
# ggplot(fin.graph, aes(stand.id, ple)) + geom_bar(stat="identity")

