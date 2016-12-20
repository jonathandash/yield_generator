#library(multicore)




### Generate N ordered block identifiers
### For ceiling(N/block.size) blocks
### So that the number in each block is approximately equal to block.size
### value = integer vector of length N
generate.block.ids <- function(N, block.size)
  {
    n.blocks <- ceiling(N / block.size) ## Number of blocks
    computed.block.size <- floor(N / n.blocks) ## Fixed number in each block
    ## Each id gets computed.block.size then the remainder are distributed randomly with replacement
    sort( c(rep(1:n.blocks,each=computed.block.size) 
            , sample(1:n.blocks, N - n.blocks * computed.block.size, repl=T)))

  }

### Get the 2d 1-based indices into the lower triangle from a 1d 1-based index
### 0  0  0  0
### 1  0  0  0
### 2  3  0  0
### 4  5  6  0
###
### e.g. get.triangular.matrix(1:6) =>
### 2,1
### 3,1
### 3,2
### 4,1
### 4,2
### 4,3
### Used to sample from triangle without replacement
get.triangular.index <- function(index)
  {
    i <- index - 1  # Convert to 0 based index
    row = floor(-0.5 + sqrt(0.25 + 2 * i))
    triangularNumber = row * (row + 1) / 2
    column = i - triangularNumber
    cbind(i=row+2,j=column+1) # Convert back to 1 based indexes in lower triangle
  }

### Take a sample of size sample.size from the lower triangle of a square matrix
### of size c(matrix.size,matrix.size)
### value = sample.size * 2 index matrix into lower triangle
### if sample.size > size of lower triangle then set of indexes representing complete lower triangle
### If matrix.size * (matrix.size - 1) / 2 is < .Machine$integer.max then
### sampling is without replacement and nrow(value) == sample.size
### Otherwise, sampling is with replacement and rows where i == j are discarded,
### leaving the sample slightly smaller than requested
sample.triangle <- function(matrix.size, sample.size)
  {
    N <- matrix.size
    triangle.size <- N * (N-1) / 2
    if ( triangle.size <= sample.size )
      {
        ## i is the row index
        ix <- subset(expand.grid(i = 1:N, j = 1:N), i > j)
      } else {
        ### If we can then sample without replacement
        if ( triangle.size < .Machine$integer.max )
          {
            ilt <- sample.int(triangle.size,sample.size,replace=F)
            ix <- get.triangular.index(ilt)
          } else {
            ### N is bigger than max int so sample row and column separately
            ### Swap sample so that i > j
            ### Discard where i == j (diagonal)
            ### This is not best solution because expected length of ix is 1-N/N^2 * sample size
            ### and this becomes sampling with replacement
            ### Fix this by switching to R V 3 and using sample.int on 64 bit ints
            ### Sampling without replacement is probably not a huge issue when the triangle size is really, really big.
            i <- sample.int(N,sample.size,replace=T)
            j <- sample.int(N,sample.size,replace=T)
            ix <- matrix( c( ifelse(i > j, i, j), ifelse(i > j, j, i)),ncol=2)
            
            ix <- ix[ i != j, ]
          }
      }
    ix
  }




### Calculate the sum of the elements within a single block of the  covariance matrix
### tki is Ni*k matrix.  Each row represents a target pixel (i). The kth column contains an index to the kth nearest reference pixel
### sigma.i is Ni*1 vector of the standard deviation of the Y value of the k nearest pixels to target pixel i
### tkj is Nj*k matrix.  Each row represents a target pixel (y). The kth column contains an index to the kth nearest reference pixel
### sigma.j is Nj*1 vector of the standard deviation of the Y value of the k nearest pixels to target pixel j
### rho.ij is an n*n matrix of the correlation between reference pixel li and reference pixel lj of variable Y
### Generates a block of the covariance matrix Ni*Nj and returns the sum of the variances and covariances within and the count within
### Values is named vector c(var,N)
knn.var.block <- function(tki, sigma.i, tkj, sigma.j, rho.ij)
  {
    k <- ncol(tki)
    ## Covariance matrix calculated as if correlation=1 for all cells
    cov.ij <- outer(sigma.i,sigma.j)

    ## Sum of correlations from 8a
    ## Generate N x k x N x k matrix of rho values
    ## rho values come from pre-calculated correlations between reference pixels
    ## Sum across 1st and third dimensions giving N*N matrix
    pij <- apply(outer(tki,tkj,FUN=function(i,j) { rho.ij[ matrix(c(i,j),ncol=2)  ]} ),c(1,3),sum)
    ## 14b. return sum of covariances within this block 
    variance <- sum(cov.ij * pij/k^2) 
    NN <- length(cov.ij)
    rm(pij, cov.ij) # tidy-up
    gc()
    
    c(var=variance,N=NN) # Return results
  }



cor.ij.fn <- function(i,j,k,tki,rho.ij,distance,rho.fn)
  {
      rho.li.lj <- outer( tki[ i,], tki[j,], FUN=function(li,lj) {rho.ij[ matrix(c(li,lj),ncol=2) ] }) 
      rho.i.lj <- rho.fn(distance[ i, tki[j,] ] ) 
      rho.li.j <- rho.fn(distance[ j, tki[i,] ] ) # 1 * k
      ( k^2 + sum(rho.li.lj) - k * (sum(rho.i.lj) + sum(rho.li.j)))
  }

### tki is N*k matrix.  Each row represents a target pixel (i). The kth column contains an index to the kth nearest reference pixel
### sigma is N*1 vector of the standard deviation of the Y value of the k nearest pixels to target pixel i
### rho.ij is an n*n matrix of the correlation between reference pixel li and reference pixel lj of variable Y
### ij => ns*2 index representing a sample of the N*N covariance matrix
### distance is N*n  matrix of geometric distances between N target pixels and n reference plots
### rho.fn is a function that calculates correlation given distance: correlation = rho.fn(distance)
knn.var.index <- function(tki, sigma, rho.ij, ij, distance, rho.fn)
  {
    k <- ncol(tki)

    ## 14b. return sum of covariances for sample of cells in covariance matrix
    ## For each index in sample
    ## Get sum of k*k correlations from rho.ij and multiply by sigma[i]*sigma[j]
    ## rho.ij values come from pre-calculated correlations between reference pixels
    ## cov is ns*1
    if ( is.null(distance) || is.null(rho.fn) )
      {
    cov.ij <- apply(ij, 1, FUN=function(x) {
      sigma[x[1]] * sigma[x[2]] *
        sum( outer( tki[ x[1],], tki[x[2],], FUN=function(li,lj) {
        rho.ij[ matrix(c(li,lj),ncol=2) ] })) 
    })
  } else {
    stopifnot(class(rho.fn) == 'function')
    ## 15b Var(YM2)
    cov.ij <- apply(ij, 1, FUN=function(x) {
      sigma[x[1]] * sigma[x[2]] * cor.ij.fn(x[1],x[2],k,tki,rho.ij,distance,rho.fn)
    })
  }
    sum(cov.ij)/k^2
    
  }


### Calculate variance of mean
### This process, unless divided, must generate an N*N covariance matrix
### with each cell calculated using a k*k correlation matrix
### So the full matrix size is potentially N*N*k*k which is stupendously big
### So this function divides the covariance matrix into manageable blocks, calculates block sums
### then discards the block
### McRoberts et al used sampling within the covariance matrix to achieve the same thing
### but I am currently working on a much smaller dataset so can process the whole covariance matrix
### in reasonable time, if I tackle it in chunks
### Each chunk represents two blocks of target pixels
### For diagonal blocks of the covariance matrix these are the same blocks of target pixels
### On the positive side, dividing the covariance matrix allows the work to be spread over multiple CPU cores
### tki is N*k matrix with each row representing a target pixel and each column representing the kth nearest reference pixel
### y is a vector of n observed y values, one for each reference pixel, ordered by reference pixel id
### rho.ij is a n x n matrix of pre-calculated spatial correlations between reference pixels for value y
### rho.ij[ i, j] is the correlation between the ith reference pixel and the jth reference pixel w.r.t. y
### method is how to process the covariance matrix
###     'block' => Process entire covariance matrix by splitting into blocks
###     'sample' => Sample the covariance matrix.  The diagonal is processed then a sample of size sample.size is taken from the lower triangle
### block.size size blocks that covariance matrix is split into for processing. This is size of 1 dimension so 500 => 500*500 cells
### sample.size => the maximum number of cells to sample from the lower triangle. If sample.size > N*(N-1)/2 then the whole matrix is processed (without blocking)
### distance is N*n  matrix of geometric distances between N target pixels and n reference plots
### rho.fn is a function that calculates correlation given distance: correlation = rho.fn(distance)
### Value is the variance of y calculated using McRoberts et al, 2007, equation 14b
knn.var <- function(tki, y, rho.ij,method,block.size,sample.size,use.mc,distance,rho.fn)
  {
    stopifnot( method  %in% c('block','sample'))
    k <- ncol(tki) # Number of nearest neighbours
    N <- nrow(tki) # Number of target pixels

    ## Nxk matrix of Y values
    ## with one row per target pixel and one col per reference pixel
    ## Cell values contain Y for that containing Y
    ## tm1 has same dimensions as tki
    tky <- matrix( y[ as.vector(tki) ], nrow=nrow(tki), ncol=ncol(tki) ) 
    if ( nrow(tky) <= 1 )
      {
        v <- NaN
      } else  { ## Calcuate variance
        ## diag.pij is    
        ## N*1 vector of the sum of correlations between reference pixels for a single target pixel
        ## diag.pij is the Sum(Sum(rho_j1,j2)) in the denominator from equation 6a
        ## And is equal to the diagonal of the correlation summation in 8a
        ## For each row in tki, create a 2d index using all combinations of the 1d reference pixel ids in that row, use this to extract the k*k correlations from rho.ij and sum these 
        diag.pij <- apply(tki, 1 , FUN=function(x) { sum(outer(x,x,FUN=function(i,j) {rho.ij[ matrix(c(i,j),ncol=2)  ]} )) } )
        
        ## N*1 vector of sigma at each target pixel based on k nearest neighbours
        ## Equation 6a
        sigma <- sqrt(apply(tky, 1, function(x) { sum((x - mean(x))^2) })  / ( k - diag.pij/k))

        if (method == 'block')
          {
            block.id <- generate.block.ids(N,block.size) # N*1 grouping variable
            n.blocks <- length(unique(block.id)) # Number of blocks
            
            ## Pairs of block ids
            pairs <- expand.grid(i=1:n.blocks,j=1:n.blocks)
            ## Because of symmetry we only need to do off-diagonal blocks once
            pairs <- with(pairs, pairs[ i <= j, ])
            pairs$weight <- with(pairs,ifelse(i == j, 1.0, 2.0)) # Weight diagonal blocks by 1 and off-diagonal by 2

            ## Select version of lapply.  If using multiple cores then use mclapply.  This won't work under Windows
            if ( use.mc ) {flapply <- mclapply} else {flapply <- lapply}
            ## For each block construct a part of the covariance matrix and sum its contents
            ## Variance is sum across all blocks
            ## block.results is pn*2 matrix, one row for each pair of blocks. First column is sum of covariances, second is number of cells
            block.results <- t(matrix(unlist(flapply(1:nrow(pairs), function(bi) {
              s1 <- block.id == pairs[bi,1]
              s2 <- block.id == pairs[bi,2]
              knn.var.block( tki[s1,], sigma[s1], tki[s2,], sigma[s2], rho.ij )
            })),nrow=2))
            w <- matrix(rep(pairs$weight,2),ncol=2) # Duplicate weights because there are two columns
            tmp <- colSums( block.results * w )
            v <- tmp[1]/tmp[2]  ### Sum of covariances / sum of block sizes
          } else { # method == 'sample'
            ii <- matrix(rep(1:N,2),ncol=2) ## matrix index for diagonal
            vii <- knn.var.index( tki, sigma, rho.ij, ii, distance, rho.fn )
            ij <- sample.triangle(N, sample.size)
            expansion <-   ( N*(N-1) )/ nrow(ij) ### off-diagonal size / sample size.  Use nrow(ij) in case sample.size was bigger than size of lower.triangle
            
            vij <- knn.var.index( tki, sigma, rho.ij, ij, distance, rho.fn ) * expansion
            v <- ( vii + vij ) / N^2 ## N^2 = N * (N-1) + N
          }
      }
    v
    
  }



### tki is N*k matrix with each row representing a target pixel and each column representing the kth nearest reference pixel
### y is a vector of observed y values, one for each reference pixel, ordered by reference pixel id and indexed by reference pixel id
### rho.ij is a n x n matrix of pre-calculated spatial correlations between reference pixels for value y
### rho.ij[ i, j] is the correlation between the ith reference pixel and the jth reference pixel w.r.t. y
### method how to process the covariance matrix
###     'complete' => Process entire covariance matrix by splitting into blocks
###     'sample' => Sample the covariance matrix.  The diagonal is processed then a sample of size sample.size is taken from the lower triangle
### block.size size blocks that covariance matrix is split into for processing. This is size of 1 dimension so 500 => 500*500 cells
### sample.size => the maximum number of cells to sample from the lower triangle. If sample.size > N*(N-1)/2 then the whole matrix is processed (without blocking)
### distance is N*n  matrix of geometric distances between N target pixels and n reference plots.  Only used for McRoberts method 2 (equation 15b).  In most cases you don't want to use this so set to NULL.
### rho.fn is a function that calculates correlation given distance: correlation = rho.fn(distance). Only used for McRoberts method 2 (equation 15b). 
### values in tki are numeric indices into both y and rho.ij
### value is a named vector containing N,k,n,nr,mu,var
knn <- function( tki, y, rho.ij, method='complete', block.size=500, use.mc=T, sample.size = 500*499/2,distance=NULL,rho.fn=NULL  )
  {
    k <- ncol(tki)


    ## Nxk matrix of Y values
    ## with one row per target pixel and one col per reference pixel
    ## Cell values contain values of Y from the reference pixel
    tky <- matrix( y[ as.vector(tki) ], nrow=nrow(tki), ncol=ncol(tki) ) 

    ## N*1 vector of mean of k nearest neighbour values for each target pixel
    mu.i <- apply(tky, 1, mean)


    ## Return mean and variance of y for all N target pixels
    v <- NA
    tmp <- try(knn.var(tki,y,rho.ij,method,block.size,sample.size,use.mc,distance,rho.fn), silent=T)
    if (class(tmp) == "numeric") v <- tmp
    c(N=nrow(tki), k=k, n=length(y), nr=length(unique(as.vector(tki))), mu=mean(mu.i), var=v)
    
  }

