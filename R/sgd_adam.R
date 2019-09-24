#' SGD-Adam
#' 
#' @param delta  event indicator
#' @param z Covariate matrix
#' @param time observed event time
#' @param knot number of basis functions for time-varing effects
#' @param facility strata
#' @param M_stop maximum stopping iterations
#' @param tol the convergence threshlod
#' @param rate
#' 
#' @return return a list of the following values
#' *theta : 
#' *time :
#' *test_SGD_all :
#' @export
#' 
#' 



SGD  <- function(delta, z, time, knot = 10 ,facility = NULL, M_stop = 10000, tol= 10^(-6),rate=0.001  ) {
  track <- 5                       # used for later likelihood increment comparison for more stable results.      
  nfold <- 5 
  K     <- knot
  N     <- nrow(z)
  p     <- ncol(z)
  theta.SGD.all = NULL
  SGD.time =  NULL
  test_SGD_all = NULL
  theta_SGD_ADAM_all    <- NULL
  
  
  
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
  
  b_spline <- bs8
  
  data3          <- cbind(time,delta,facility, b_spline,z)
  data3          <- data3[order(facility),]
  m        <- as.vector(sapply(split(data3[,3],factor(data3[,3])),length))  #record the length of each strata, used for split the data.
  
  
  m_stratify     <- 0
  theta_stratify <- matrix( rep(0, knot*p), nrow=p)  #dim P*knot
  
  
  likelihood_SGD_ADAM_all <- dloglik_likelihood_gradient(knot,facility,delta,z,b_spline,theta_stratify)   #calculate the likelihood
  G              <- rep(0,p*K)
  M              <- rep(0,p*K)
  gamma=0.01 #0.01 performs better (better than 1, 0.1, 0.001, 0.0001), and is recommended by https://arxiv.org/pdf/1609.04747.pdf
  
  
  #start the loop
  stage2_key=FALSE
  while (stage2_key==FALSE){
    m_stratify <- m_stratify+1
    old_time   <- proc.time()
    
    #splite the data, select mini-batch to update theta. (proportionly from each strata.)
    bootsample3<- stra_sampling(m,nfold)
    data_sub4  <- data3[bootsample3,]
    
    time_sub4  <- data_sub4[,1]
    delta_sub4 <- data_sub4[,2]
    facility_sub4 <- data_sub4[,3]
    b_spline_sub4 <- data_sub4[,4:(knot+3)]
    z_sub4        <- data_sub4[,-(1:(knot+3))]
    
    delta_sub4    <- delta_sub4[order(time_sub4)]
    facility_sub4 <- facility_sub4[order(time_sub4)]
    z_sub4        <- z_sub4[order(time_sub4),]
    b_spline_sub4 <- b_spline_sub4[order(time_sub4),]
    time_sub4     <- time_sub4[order(time_sub4)]
    n_sub4        <- length(facility_sub4)
    
    #calculate the first order derivative. This function could do more, need further ...
    result=GDboost_gradient(knot,rate,facility_sub4,delta_sub4,z_sub4,b_spline_sub4,theta_stratify)
    
    G_temp=rep(0,p*K)
    M_temp=rep(0,p*K)
    M_temp=result$L1
    G_temp=result$L1**2
    
    G <- 0.999*G + 0.001*G_temp
    M <- 0.9*M + 0.1*M_temp
    
    M_hat <- M/(1-0.9^m_stratify)     
    G_hat <- G/(1-0.999^m_stratify)
    
    #update theta
    if (m_stratify>1){
      theta_stratify <- theta_stratify-matrix(gamma/(sqrt(G_hat)+10^-8)*M_hat,nrow=p,byrow = TRUE)
    }
    if (m_stratify==1){
      theta_stratify <- theta_stratify-gamma*matrix(result$L1,nrow=p,byrow = TRUE)
    }
    theta_SGD_ADAM_all[[m_stratify]] <- theta_stratify
    
    
    likelihood_stratify     <- dloglik_likelihood_stratify(knot,facility,delta,z,b_spline,theta_stratify)  #calculate the likelihood.
    likelihood_SGD_ADAM_all <- c(likelihood_SGD_ADAM_all, likelihood_stratify)
    

    # if the relative likelihood change is small within certain time, we break the loop
    if (m_stratify>=(10+track)) {
      llk.diff.all=NULL
      for (key in 1:track) {
        llk.diff.all = c(llk.diff.all, likelihood_SGD_ADAM_all[m_stratify-key+1]-likelihood_SGD_ADAM_all[m_stratify-key])
      }
      llk.diff=max(llk.diff.all)
      llk.diff_null = likelihood_SGD_ADAM_all[m_stratify]-likelihood_SGD_ADAM_all[1]
      dist=theta_SGD_ADAM_all[[m_stratify]]-theta_SGD_ADAM_all[[m_stratify-1]]
      if(abs(llk.diff/llk.diff_null) < tol) {
        stage2_key=TRUE 
        break
      }
    }
    
    if (m_stratify==M_stop){
      stage2_key=TRUE
      break
    }
  } #end while
  
  
  theta.SGD.all <- rbind(theta.SGD.all, theta_stratify)
  
  temp          <- ddloglik(knot,facility,delta,z,bs8,theta_stratify,length(unique(facility)))
  
  test_SGD      <- rep(0,p)
  constrain     <- -diff(diag(knot*1),differences=1)
  j=0
  repeat{
    j=j+1
    
    theta_SGD_j= theta_stratify[j,]
    
    L2=solve(temp$GVG)[((j-1)*knot+1):((j-1)*knot+knot),((j-1)*knot+1):((j-1)*knot+knot)]
    test_contrast=t(constrain%*%theta_SGD_j)%*%solve(constrain%*%L2%*%t(constrain))%*%(constrain%*%theta_SGD_j)
    test_SGD[j]=1-pchisq(test_contrast, (knot-1))
    if(j==p) break
  }
  
  test_SGD_all  <- rbind(test_SGD_all, test_SGD)
  
  return( list(theta = theta.SGD.all, time = SGD.time, test_SGD_all = test_SGD_all ))
}













stra_sampling<-function(size, prop){
  m     <- length(size)
  mark  <- cumsum(size)
  total <- sum(size)
  sample.id <- NULL
  sample.id <- c(sample.id, sample.int(size[1], size=floor(size[1]/prop), replace=FALSE))
  if (m>1) {
    for(i in 2:m){
      sample.id <- c(sample.id,mark[i-1]+sample.int(size[i],size=floor(size[i]/prop),replace=FALSE))
    }
  }
  return(sample.id)
}


strata <- function(n,nfolds){
  facility=rep(0,n)
  temp <- split(sample(1:n), rep(1:nfolds, length = n))
  for (i in 1:nfolds){
    facility[temp[[i]]]=i
  } 
  return(facility)
}








