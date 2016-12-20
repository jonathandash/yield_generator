library(gstat)


### Functions for generating a correlation matrix and local variance from the correlation matrix

### Extract empirical semivariogram function from variogramModel
### Function calculates gamma(d) where d is distance and gamma is the semi-variogram
### Only works where model = 'Exp', with or without nugget
exp2f <- function(vf)
  {
    stopifnot( vf$model[2] == 'Exp' && vf$range[1] == 0)
    function(x)
      {
        vf$psill[1] + vf$psill[2] * ( 1 - exp(-x / vf$range[2] ))
      }
    
  }


### Calculate average and variance of k nearest neighbours to a point
### rho.ij is n*n correlation matrix between the k-nearest neighbours of the reference points, ordered by reference point identifier.
### yi is n*1 vector of observed y values for each reference point, ordered by reference point identifier
### idk is n*k matrix of the reference point identifiers for the k nearest neighbours for each of n reference points, ordered by reference point identifier
### Caution:
### It is assumed that reference point identifiers are complete 1:n and the identifiers in idk may be used as integer indexes into idk, yi and rho.ij (both dimensions).
### value = list(yknn.i, sigma.i)
### yknn.i is n*1 vector where each element is the unweighted average of the k values from yi indexed by the k indexes in each row of idk  
### sigma.i is n*1 calculated using equation 6a of McRoberts et al, 2007
yknn.i <- function(rho.ij, yi, idk)
  {
    k <- ncol(idk)
    yik <- matrix( yi[as.vector(as.matrix(idk))], ncol=k) # n*k matrix of y values 
    rownames(yik) <- rownames(yi)
    yknn <- rowMeans(yik) ### Knn estimate for reference point equal to mean of reference values

    ## n*1 Sums of correlations for the k-nearest neighbours for each of the n reference points
    ## For each row, create k*k matrix containing correlation between each pair of reference points, then sum
    diag.pij <- apply(idk, 1 , FUN=function(x) { sum(outer(x,x,FUN=function(i,j) {rho.ij[ matrix(c(i,j),ncol=2)  ]} )) } )
    var.yknn <- apply(  (yik -  yknn)^2, 1, sum) / ( k - diag.pij/k ) # Equation 6a, McRoberts 2007
    list(yknn.i=yknn, sigma.i=sqrt(var.yknn))
  }

### Implements one iteration of equations 6a-7b of McRoberts et al, 2007
### n is number of reference points.
### k is number of nearest neighbours
### rho.ij is an initial or intermediate n*n correlation matrix between the k-nearest neighbours of the reference points, ordered by reference point identifier in both dimensions
### yi is n*1 vector of observed y values for each reference point, ordered by reference point identifier
### idk is n*k matrix of the reference point identifiers for the k nearest neighbours for each of n reference points, ordered by reference point identifier
### xy is n*2 matrix (or data.frame) containing x,y co-ordinates of n reference pixels, ordered by reference pixel identifier
### Caution:
### It is assumed that reference point identifiers are complete 1:n and the identifiers in idk may be used as integer indexes into xy, idk, yi and rho.ij (both dimensions).
### cutoff is limit above which distances are not used to build variogram (see gstat::variogram)
### value is list containing new version of rho.ij and var.i, the variance at each of the n reference points of the k nearest neighbours
### Uses unweighted means and variances
improve.correlation.matrix <- function(rho.ij, yi, idk, xy,cutoff)
  {
    k <- ncol(idk)
    yik <- matrix( yi[ as.vector(as.matrix(idk)) ], ncol=k)
    tmp <- yknn.i(rho.ij, yi, idk)
    yknn <- tmp$yknn.i; sigma.yknn <- tmp$sigma.i
    delta <- (as.vector(yi) - yknn) / sigma.yknn
    s <- ! is.na(delta) 
    rho.ij.2 <- rho.ij
    vf <- NULL; gamma.fn <- NULL;
    vd <- try(variogram(delta[s]~1,locations=formula(paste('~',paste(names(xy),collapse="+"),sep='')),xy[s,],cressie=T,cutoff=cutoff),silent=T)
    if ("gstatVariogram" %in% class(vd))
      {
        # Use na.rm=T because delta == NaN if sigma == 0
        vf <- try(fit.variogram(vd,model=vgm(max(delta,na.rm=T),"Exp",cutoff,max(delta,na.rm=T)/2)),silent=T)
        if ( ! is.null(vf) && "variogramModel" %in% class(vf) &&  ! is.null(attr(vf,'singular')) && ! attr(vf,'singular') )
          {
            ## Generate another correlation matrix 
            gamma.fn <- exp2f(vf)  # Eq 7a semi-variogram function
            dij <- as.matrix(dist(xy)) # dij for equation 7b  n*n distance matrix
            infinity <- max(dij) * 2   # gamm atotal from 7b
            rho.ij.2 <- 1 - gamma.fn(dij) / gamma.fn(infinity) # Equation 7b
            diag(rho.ij.2) <- 1 # Set diagonals to 1 because if there is a nugget, they won't be 1
          }
      }
    list(rho.ij=rho.ij.2, vd=vd,vf=vf,gamma.fn=gamma.fn)
    
  }


### Implements multiple iteration of equations 6a-7b of McRoberts et al, 2007
### n is number of reference points.
### k is number of nearest neighbours
### yi is n*1  vector of observed y values for each reference point, ordered by reference point identifier
### idk is n*k matrix of the reference point identifiers for the k nearest neighbours for each of n reference points, ordered by reference point identifier
### xy is n*2 matrix (or data.frame) containing x,y co-ordinates of n reference pixels, ordered by reference pixel identifier
### Caution:
### It is assumed that reference point identifiers are complete 1:n and the identifiers in idk may be used as integer indexes into xy, idk and yi
### Assumes that yi, xy and idk are ordered the same; ie that each row of each input matrix represents the same reference point
### cutoff is limit above which distances are not used to build variogram (see gstat::variogram)
### tolerance is maximum acceptable difference between elements of successive approximation of the correlation matrix
### value is n*ncorrelation matrix
### Uses unweighted means and variances
calculate.correlation.matrix <- function(yi, idk, xy,cutoff=1000,tolerance=0.01,maxiter=50)
  {
    stopifnot("numeric" %in% class(yi) && length(yi) == nrow(idk) && length(yi) == nrow(xy) && ncol(xy) == 2)
    stopifnot(all(rownames(xy) == rownames(idk)))
    nr <- length(yi)
    rho.ij.new <- matrix(0,ncol=nr,nrow=nr) ## First attempt at correlation matrix assumes no spatial correlation
    diag(rho.ij.new) <- 1

    happy <- FALSE
    n <- 0
    # iterate until 
    while ( ! happy && n < maxiter )
      {
        rho.ij <- rho.ij.new
        n <- n + 1
        tmp <- improve.correlation.matrix(rho.ij, yi, idk, xy,cutoff)
        rho.ij.new <- tmp$rho.ij
        happy = all( abs(rho.ij.new - rho.ij) < tolerance)
      }
    if (n >= maxiter)
      {
        warning(paste("Maximum iterations (",maxiter,") reached with maximum difference between successive matrices equal to ",max(abs(rho.ij.new - rho.ij)),sep=""))
      }
    list(rho.ij=rho.ij.new,n=n,vd=tmp$vd,vf=tmp$vf)
    
  }
