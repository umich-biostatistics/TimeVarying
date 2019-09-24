#' Netwon Raphson
#' 
#' @param delta  event indicator.
#' @param z Covariate matrix.
#' @param time observed event time.
#' @param knot number of basis functions for time-varing effects.
#' @param facility strata.
#' @param M_stop maximum stopping iterations.
#' @param tol the convergence threshlod.
#' @param rate
#' 
#' @return return a list of the following values
#' *theta : 
#' *test_SGD_all :
#' *b_spline :
#' 
#' @examples
#' 
#' # the SGD
#' # generate the simuluation data
#' data2 <- simul(N = 2000, N_Strata = 5, p=5 )
#' result2 <- SGD(data2$delta, data2$z, data2$time, facility = data2$facility,
#'                knot = 10, M_stop = 10000, tol= 10^(-6),rate=0.001)
#' 
#' @export
#' 
#' 

NR_new <- function(delta, z, time, knot = 10 ,facility = NULL, M_stop = 10000, tol= 10^(-6),rate=0.001  ){
  old_time <- proc.time() # record used time
  track <- 5
  K     <- knot
  diff_0<- 0
  MSE_0 <- 0
  N     <- nrow(z)
  p     <- ncol(z)
  
  NR.time      <- NULL    # record time used
  theta.NR.all <- NULL    # record theta
  test_NR_all  <- NULL    # record test result
  
  #order the variables
  delta    <- delta[order(time)]
  facility <- facility[order(time)]
  z        <- z[order(time),]
  time     <- time[order(time)]
  
  #using bspline to span the function.
  time2    <- time[delta==1]      # record the event time. 
  knot_set <- quantile(time2,prob=seq(1:(knot-4))/(knot-3)) 
  bs7      <- splines::bs(time,df=knot, knot=knot_set, intercept=TRUE, degree=3)
  bs8      <- matrix(bs7, nrow=N) # for build beta_t
  
  #record theta
  theta_NR <- matrix(rep(rep(0,p), knot), nrow=p, byrow=FALSE)
  
  number_facility = length(unique(facility))
  nrloop   <- NRloop(knot, facility, delta ,z, bs8, theta_NR, M_stop, number_facility, rate, tol)
  
  theta_NR <- nrloop$theta
  likelihood_NR_all <- nrloop$likelihood_NR_all

  
  return(list(theta = theta_NR,
              likelihood = likelihood_NR_all,
              b_spline    = bs8))
  

  
}
