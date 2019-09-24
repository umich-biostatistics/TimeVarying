
#' This simulation function is used to simulate data for testing 
#' the functions in this package.
#' 
#' @param N the sample size for each strata.
#' @param N_Strata number of stratum.
#' @param p number of parameters.
#' @param p_true number of true parameters.
#' 
#' @return return a list of the following values
#' *delta :  event indicator
#' *z: Covariate matrix
#' *facility: 
#' *time: the death time
#' 
#' @examples
#' # generate the simuluation data
#' data2 <- simul(N = 2000, N_Strata = 5, p=5 )
#' 
#' @export
#' 

simul <- function(N = 1000, N_Strata = 10, p=5 ){
  
  n_f   <-  rpois(N_Strata, lambda = N)      
  N     <-  sum(n_f)                         # total number of samples
  
  gamma <- rnorm(N_Strata, mean=0, sd=0.5)   
  gamma_subject <- rep(gamma,n_f)

  F_pre    <- 1:N_Strata
  facility <- rep(F_pre, n_f)                # record the strata
  
  Sigma_z1 <- AR1(0.6,p)                     # compute the sigma
  
  z        <- rmvnorm(N, mean=rep(0,p), sigma=Sigma_z1) 
  z_012_rare=function(x){
    # U=runif(2, 0.97,1)
    U=runif(2, 0.85,0.95)
    U=U[order(U)]
    
    x2=quantile(x,prob=U)
    x3=x
    x3[x<x2[1]]=0
    x3[x>x2[2]]=2
    x_index=(((x<x2[1])+(x>x2[2]))<1)
    x3[x_index]=1
    return(x3)
  }
  
  z    <- apply(z,2, z_012_rare)
  
  MAF  <- apply(z,2,mean)
  
  U    <- runif(N, 0,1)
  
  pre_time <- rep(0, N)
  for (i in 1:(N)) {
    f=function(t) {
      #integrand <- function(x) {0.5*exp(gamma_subject[i]+z[i,1]-z[i,3]+sin(3*pi*x/4)*(x<3)*z[i,2]+exp(-x)*(x<3)*z[i,4]+z[i,5])}
      integrand <- function(x) {0.5*exp(gamma_subject[i]+z[i,1]-z[i,3]+sin(3*pi*x/4)*(x<3)*z[i,2]-z[i,4]+z[i,5])}
      
      Lambda=integrate(integrand, lower = 0, upper = t)$value
      Lambda+log(1-U[i])
    }
    r1 <- suppressWarnings(try(uniroot(f,  lower = 0, upper = 4), silent=TRUE))
    if (class(r1) == "try-error"){    
      pre_time[i]=4
    }
    else pre_time[i]=uniroot(f,  lower = 0, upper = 4)$root
  }
  
  pre_censoring  <- runif(N,0,3)
  pre_censoring  <- pre_censoring*(pre_censoring<3) + 3*(pre_censoring>=3)
  tcens <- (pre_censoring<pre_time) # censoring indicator
  delta <- 1-tcens
  time  <- pre_time*(delta==1) + pre_censoring*(delta==0)
  
  return(list( delta=delta, z=z, facility = facility, time=time))
  
}


AR1 <- function(tau, m) {
  if(m==1) {R <- 1}
  if(m > 1) {
    R <- diag(1, m)
    for(i in 1:(m-1)) {
      for(j in (i+1):m) {
        R[i,j] <- R[j,i] <- tau^(abs(i-j))
      }
    }
  }
  return(R)
}



