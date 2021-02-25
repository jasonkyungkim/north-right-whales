# Functions




solow <- function(t, pi=0.5, now)  {
  # t: a vector of sighting dates (assumed in decimal years)
  # return: probability that an individual is still alive
  # pi: prior probability
  input <- t
  end_time <- now
  
  t <- sort(t)
  n <- length(t) - 1
  T <- now - t[1]
  t <- t - t[1]
  B <- (n-1)/((T/t[n+1])^(n-1)-1)
  # three or more sightings
  if(n>1)  return ( 1/(1+((1-pi)/(pi*B))))

  #two sightings or less
  # Caley function
  #if(n<=1)  {
  #  whale_sightings <- convert_input(input, month=30.4375, start_time=0 ,end_time=end_time)
  #  nonconstant <- fit.func.nonconstant(iter=n.iter,          # non-constant parameters
  #                                     y=whale_sightings,    #sighting data ex. (0,1,1,0,0)
  #                                     pgr.init=0.0,         #population growth rate
  #                                     delta.init=0.69,      #yearly detection rate per unit
  #                                     eps0.init=0,          #epsilon 1
  #                                     eps1.init=1,          #epsilon 2
  #                                     N0.init=1)            #population size
  # nonconstant <- nonconstant[-(1:burn.in),]
  # return (mean(is.na(nonconstant[,6])))
  #}
  
   if(n==1)  {
     gap <- t[2]-t[1]
     if(T<gap)      prob <- .7
     if(T<gap-3)    prob <- .9
     if(T==gap)     prob <- .5
     if(T>gap)      prob <- .3
     if(T>gap+3)    prob <- .1
     return(prob)
   }
   # one sighting
   if(n==0)  {
     if(T<=1)  prob <- .9
     if(T==2)  prob <- .7
     if(T==3)  prob <- .5
     if(T==4)  prob <- .3
     if(T==5)  prob <- .1
     if(T>=6)  prob <- 0
     return(prob)
   }
}


converttime <- function(year, month, day, time)  {
  # convert database time format to standard time format
  minute <- time%%100
  hour <- (time-minute)/100
  sightingTime <- ISOdate(year, month, day, hour, minute)
  base <- ISOdate(0, 1, 1, 00, 00, 00)
  return(difftime(sightingTime, base))
}




# calculates the geodesic distance between two points specified by radian lat/long using the
#   Spherical Law of Cosines (slc)
# adapted from http://www.r-bloggers.com/great-circle-distance-calculations-in-r/
gcd.slc <- function(long1, lat1, long2, lat2) {
  if(identical(lat1,lat2) & identical(long1,long2))  {
    return(0)
  }  else  {
    R <- 6371            # Earth mean radius [km]
    d <- acos(sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2) * cos(long2-long1)) * R
    return(d)            # distance in km
  }
}



# convert degrees to radians
deg2rad <- function(deg) return(deg*pi/180)



invlogit <- function(x)  {
  return( 1 / (1 + exp(-x)) )
}



logit <- function(p)  {
  return( log(p/(1-p)) )
}

#month=30.49657
convert_input <- function(input, month=30.4375, start_time=0, end_time=now) 
{
  #Function that converts sighting dates into a vector of 0 (not seen) and 1's (seen)
  #for every month starting from the first sighting date to the current time
  #Parameters
  # input - vector of a whale's sighting records
  # month - the duration of a month in the expected time format (julian)
  # start_time - number of months before the initial sighting to begin the sighting records
  # end_time - Current time stamp in the expected time format
  #Return
  # input_sightings - vector for every month before initial sighting start and current time stamp
  #                   with 0's or 1's representing not seen or seen for that whale.
  
  input_sightings <- rep(0, start_time)
  sighting = 1
  curr_time <- input[sighting]                   # initializing time to first sighting date
  done = FALSE
  
  # Realistically should be reversing the order that the time is being incremented, starting with now variable time and decrementing down to the first sighting
  while (!done) { 
    if (input[sighting] <= curr_time & sighting <= length(input)) {     # Assuming we already check that sightings are all independent
      input_sightings <- append(input_sightings, 1)       # Whale was seen, add 1 to vector
      sighting <- sighting + 1                            # increment sighting index of focus
    } else {
      input_sightings <- append(input_sightings, 0)       # Whale was not seen, add 0 to vector
    }
    if (curr_time > now) done = TRUE
    
    curr_time <- curr_time + month               # incrementing the current month where we are checking if the whale was seen or not
  }
  input_sightings <- append(input_sightings, 0) #added in to bypass error of having last sighitng also be the last value in vector, error at line 400
  
  input_sightings #Return value
}


# Non-Constant
library(MASS)


# Set burn in and length of chain
burn.in <- 1E3
n.iter <- 1E4 + burn.in


sim.N <- function(pgr,time,N0) {
  # Returns:
  #	population trajectory from start of observation period
  # Args:
  #	N0 - initial population size
  #	pgr - exponential population growth rate per year
  #	time - time to end of observation period
  ans <- numeric(time)
  ans[1] <- N0	
  if(time>1) {
    for (t in 2:time){
      ans[t] <- ans[t-1]*exp(pgr)
    }
  } 
  ans
}

ilogit <- function(x)  {
  # Returns:
  #	inverts logit transformation	
  1/(1+exp(-x))
}

lam.gen <- function(delta,N) {
  # Returns:
  # 	vector of yearly detection probability based on population size
  # Args: 
  # 	delta -- per capita detection rate 
  # 	N -- population size
  1-exp(-delta*N)	
}

phi.gen <- function(eps0,eps1,N) {
  # Returns:
  # 	vector of extinction probability based on population size
  # Args: 
  #	eps0 -- intercept for logit of extinction probability
  #	eps1 -- population coefficient for logit of extinction probability
  
  lp <- eps0 - eps1*log(N)
  # Avoiding numerical overflow
  if (lp>50) {
    ans <- 1
  } else {
    ans <- ilogit(lp)
  }
  ans
}


correct.negative <- function(x) {
  # For correcting numerical errors for non-negative vectors
  x[sign(x)<0] <- 0
  x
}

impute.TE <- function(p.cease) {
  # Returns:
  #	 Estimated time to extinction (TE) within observation window, o.w. returns NA if right-censored
  # Args:
  # 	p.cease -- a vector of yearly estimates of extinction during the observation window (p.cease) 
  
  # get rid of any floating point errors
  p.cease <- correct.negative(p.cease)	
  
  # Probability of being right-censored
  p.cens <- tail(p.cease,1)
  if(sign(p.cens<0)) p.cens <- 0
  
  # Random draw for time to extinction conditional on p.cease
  TE.draw <- rmultinom(1,1,p.cease)
  
  # Where does extinction occur
  result <- (1:length(p.cease))[TE.draw==1]	
  if (result <= (length(p.cease)-1)) {
    return(result)
  } else {
    return(NA)
  }
}



####################################################################################################################################

calc.pars <- function(p, y=whale_sightings) {		 
  # Args: 
  # 	 p - c(pgr (pop. growth rate), delta, eps0, eps1, N0)
  # Function calls: 
  #	sim.N(), lam.gen(), phi.gen(), impute.TE() 
  # Returns:
  #	y.obs -- observations up until imputed time of extinction
  #	N -- vector containing population trajctory
  #	phis -- vector of extinction probabilities
  #	lambdas -- vector of detection probabilities	
  #	p.cease -- vector of probability of extinction over the zero.count period and beyond
  #	p.cens -- probability of right-censoring in extinction year [also present as last entry in p.cease]
  # 	TE -- imputed time of extinction from impute.TE()
  
  pgr <-  p[1]
  delta <- p[2]
  eps0 <- p[3] 
  eps1 <- p[4]
  N0 <- p[5]
  obs <- length(y)
  
  # positions of sightings
  pos <- (1:length(y))[y==1]		
  
  # time of last sighting [from beginning of observation period]
  final <- tail(pos,1)		
  
  # observations since last sighting
  zero.count <- obs - final		
  
  # projected population trajectory up until end of current observation period
  N.traj <- sim.N(N0=N0, pgr=pgr, time=obs)			
  
  # vector of detection probs
  lambdas <- sapply(N.traj, function(x) lam.gen(delta=delta,N=x))	
  
  # vector of extinction probs
  phis <- sapply(N.traj, function(x) phi.gen(eps0=eps0,eps1=eps1,N=x))	
  
  # P(extinction) by year following last sighting "+1" for right censoring
  p.cease <- numeric(obs+1)		
  
  # Can't have ceased prior to last sighting
  p.cease[1:final] <- 0			
  # First extinction probability following last sighting is simply phi for that time
  p.cease[final+1] <- phis[final+1]
  if(zero.count>1) {
    # for the remaining period of no sightings
    for (d in (final+2):obs) {	
      #log.ans <- sum(log(1-phis[(final+1):(d-1)])) + sum(log(1-lambdas[(final+1):(d-1)])) + log(phis[d])
      #p.cease[d] <- exp(log.ans)
      p.cease[d] <- prod(1-phis[(final+1):(d-1)])*prod(1-lambdas[(final+1):(d-1)])*phis[d] 	# potentially numerically challenging
    }
  }
  # get rid of any floating point errors
  p.cease <- correct.negative(p.cease)	
  # calculating probability of right-censoring for TE - Not extinct & not seen for zero.count
  if(sum(p.cease)>0) {
    log.ans <- sum(log(1-phis[(final+1):obs])) + sum(log(1-lambdas[(final+1):obs]))
    p.cease[obs+1] <- exp(log.ans)
    #p.cease[obs+1] <- prod(1-phis[(final+1):obs])*prod(1-lambdas[(final+1):obs])		# potentially numerically challenging
  } else {
    p.cease[obs+1] <- 1				# right-censored for sure as zero prob. mass for TE<T
  } 
  # select only observation data for period following last sighting
  p.cease <- p.cease[(final+1):(final+zero.count+1)]	
  
  # Normalizee
  p.cease <- p.cease/sum(p.cease)		
  p.cens <- tail(p.cease,1)
  
  # Calculate imputed extinction time [from beginning of observation period]
  TE <- final + impute.TE(p.cease)   	
  
  if(!is.na(TE)) {					
    # TE occurs during observation period
    phis <- head(phis,TE)		# Trim to length of TE
    lambdas <- head(lambdas,TE)	# Trim to length of TE
    y.obs <- head(y,TE)		# Trim to length of TE
    N.traj <- head(N.traj,TE) 	# Trim to length of TE
  } else {
    # If right-censored then all observations up until end of observation period are included
    y.obs <- y				
  }
  list(N=N.traj, p.cease=p.cease, p.cens=p.cens, TE=TE, y.obs=y.obs, phis=phis, lambdas=lambdas )
}

####################################################################################################################################


calc.pars.given.TE <- function(p,TE,y=y) {	# p=prop.p TE=TE.imp
  # Function to calculate the phis and the lambdas conditional on TE
  # Args: 
  #  N0 - intial population size
  #  p - proposed pars (pgr,delta,eps0)
  #  y - observed data
  # Function calls: sim.N(), phi.gen(), lam.gen()
  pgr <-  p[1]
  delta <- p[2]
  eps0 <- p[3] 
  eps1 <- p[4]
  N0 <- p[5]
  
  if(!is.na(TE)) {
    # population trajectory until extinction prior to end of observation period
    N.traj <- sim.N(N0=N0, pgr=pgr, time=TE)							
  } else {
    # population trajectory until right censoring
    N.traj <- sim.N(N0=N0, pgr=pgr, time=length(y))						
  }
  
  # detection probs until right censoring
  lambdas <- sapply(N.traj, function(x) lam.gen(delta=delta,N=x))				
  
  # extinction probs until right censoring 
  phis <- sapply(N.traj, function(x) phi.gen(eps0=eps0,eps1=eps1,N=x))			
  
  list(N=N.traj, TE=TE, phis=phis, lambdas=lambdas)
}

# Log likelihood for lamba
lnL.lam <- function(lams,ys){
  sum(log(lams[ys!=0])) + sum(log(1-lams[ys!=1]))
}

# Log likelihood for phi
lnL.phi <- function(phis,TE) {
  if(!is.na(TE)) {			
    # Extinct before end of observation period
    ans <- sum(log(1-phis[-length(phis)])) + log(tail(phis,1))
  } else {				
    # Right-censored (extinction beyond observation period)
    ans <- sum(log(1-phis))
  }
  ans
}

#####################################################################################################################################################

# Specify prior functions	
N0.prior <- function(x) dunif(x,5,50,log=F)				
lnL.N0.prior <- function(x) dunif(x,5,50,log=T)		

pgr.prior <- function(x) dunif(x,min=-2.3, max=0.69,log=F)		 
lnL.pgr.prior <- function(x) dunif(x,min=-2.3, max=0.69,log=T)	 

delta.prior <- function(x) dunif(x,0.01,4.6,log=F)
lnL.delta.prior <- function(x) dunif(x,0.01,4.6,log=T)	 

eps0.bounds <- c(-20,20)
eps0.prior <- function(x) dunif(x,min(eps0.bounds),max(eps0.bounds),log=F)			
lnL.eps0.prior <- function(x) dunif(x,min(eps0.bounds),max(eps0.bounds),log=T)			

eps1.bounds <- c(0,20) 
eps1.prior <- function(x) dunif(x,min(eps1.bounds),max(eps1.bounds),log=F)
lnL.eps1.prior <- function(x) dunif(x,min(eps1.bounds),max(eps1.bounds),log=T)

#####################################################################################################################################################

# Proposal functions
# Proposals shared amongst models
propose.N0 <- function(x){
  return(x)}
propose.pgr <- function(x, sd=0.05) {		
  rnorm(1,x,sd=sd)}
propose.delta <- function(x,sd=0.5) {
  rlnorm(1,log(x),sd)}
q.delta <- function(x1,x2) {dlnorm(x1,log(x2),sd=0.5,log=T)} 

propose.eps0 <- function(x, sd=1.5) {
  rnorm(1,x,sd=sd)}
propose.eps1 <- function(x, sd=1.5) {
  rnorm(1,x,sd=sd)} 

propose.eps1 <- function(x,sd=0.15) {
  rlnorm(1,log(x),sd)}
q.eps1 <- function(x1,x2) {dlnorm(x1,log(x2),sd=0.15,log=T)} 

#####################################################################################################################################################

# Function to run MCMC sampler with MH 
fit.func.nonconstant <- function(N0.init=1,y=whale_sightings,iter=100, pgr.init=0.0, delta.init=NA, eps0.init=0.0, eps1.init=-0.1) {	
  
  # NB "prop" is shorthand for proposed
  # NB "curr" is shorthand for current
  
  # Set delta.init such that initial lambda is 0.5
  delta.init <- -log(1-0.5)/N0.init
  
  # Initialize current parameters
  curr.p <- c(pgr.init,delta.init,eps0.init,eps1.init,N0.init)  
  
  
  # #What
  # cat(y)
  # obs <- length(y)
  # cat(obs)
  
  
  # Setup
  result <- matrix(0,nrow=iter,ncol=10)	# pgr, delta, eps0, eps1, N0, TE, acceptance indicators (x4)
  for (i in 1:iter) {	#i=1	i=i+1
    #cat("Doing interation", i, "of", iter, "\n") 
    
    # Generate values of y's & T_E by imputation conditional on current N0, delta, eps0 and eps1 
    curr.vals <- calc.pars(p=curr.p,y=y)
    
    # imputed value for TE -- stays fixed for remainder of conditional sampling
    TE.imp <- curr.vals$TE			
    
    # observed observations corresponding to TE.imp
    y.obs <- curr.vals$y.obs			
    
    # Update delta
    
    curr.lnLik.delta <- lnL.lam(lams=curr.vals$lambdas,ys=y.obs) + lnL.delta.prior(curr.p[2])
    prop.p <- curr.p
    prop.p[2] <- propose.delta(curr.p[2])
    
    prop.vals <- calc.pars.given.TE(p=prop.p,TE=TE.imp,y=y)	
    
    
    prop.lnLik.delta <- lnL.lam(lams=prop.vals$lambdas,ys=y.obs) + lnL.delta.prior(prop.p[2])
    
    LR <- exp(prop.lnLik.delta - curr.lnLik.delta + q.delta(curr.p[2],prop.p[2]) - q.delta(prop.p[2],curr.p[2]))
    
    if(runif(1) < LR & !is.na(LR)) {
      curr.p <- prop.p
      curr.vals <- prop.vals
      result[i,8] <- 1 	# indicate acceptance
    } else {
      result[i,8] <- 0 	# indicate rejection
    }
    
    # Update epsilons
    curr.lnLik.eps <-  lnL.phi(phis=curr.vals$phis,TE=TE.imp) + lnL.eps0.prior(curr.p[3]) + lnL.eps1.prior(curr.p[4]) 
    
    # Update eps0 annd eps1 jointly
    prop.p <- curr.p
    
    prop.p[3:4] <- c(propose.eps0(curr.p[3]),propose.eps1(curr.p[4]))
    
    prop.vals <- calc.pars.given.TE(p=prop.p,TE=TE.imp,y=y)	
    
    prop.lnLik.eps <- lnL.phi(phis=prop.vals$phis,TE=TE.imp) + lnL.eps0.prior(prop.p[3]) + lnL.eps1.prior(prop.p[4]) 
    
    LR <- exp(prop.lnLik.eps - curr.lnLik.eps + q.eps1(curr.p[4],prop.p[4]) - q.eps1(prop.p[4],curr.p[4]))
    
    if(runif(1) < LR & !is.na(LR)) {
      curr.p <- prop.p
      curr.vals <- prop.vals
      result[i,9] <- 1 	# indicate acceptance
    } else {
      result[i,9] <- 0 	# indicate rejection
    }
    
    # Update pgr
    curr.lnLik.pgr <- lnL.lam(lams=curr.vals$lambdas,ys=y.obs) + lnL.phi(phis=curr.vals$phis,TE=TE.imp) + lnL.pgr.prior(curr.p[1])
    prop.p <- curr.p
    prop.p[1] <- propose.pgr(curr.p[1])
    
    prop.vals <- calc.pars.given.TE(p=prop.p,TE=TE.imp,y=y)	
    
    prop.lnLik.pgr <- lnL.lam(lams=prop.vals$lambdas,ys=y.obs) + lnL.phi(phis=prop.vals$phis,TE=TE.imp) + lnL.pgr.prior(prop.p[1])
    LR <- exp(prop.lnLik.pgr - curr.lnLik.pgr)
    
    if(runif(1) < LR & !is.na(LR)) {
      curr.p <- prop.p
      curr.vals <- prop.vals
      result[i,7] <- 1 	# indicate acceptance
    } else {
      result[i,7] <- 0 	# indicate rejection
    }
    
    # Update N0
    curr.lnLik.N0 <- lnL.lam(lams=curr.vals$lambdas,ys=y.obs) + lnL.phi(phis=curr.vals$phis,TE=TE.imp) + lnL.N0.prior(curr.p[5])
    prop.p <- curr.p
    prop.p[5] <- propose.N0(curr.p[5])
    
    prop.vals <- calc.pars.given.TE(p=prop.p,TE=TE.imp,y=y)	
    
    prop.lnLik.N0 <- lnL.lam(lams=prop.vals$lambdas,ys=y.obs) + lnL.phi(phis=prop.vals$phis,TE=TE.imp) + lnL.N0.prior(prop.p[5])
    LR <- exp(prop.lnLik.N0 - curr.lnLik.N0) 
    
    if(runif(1) < LR & !is.na(LR)) {
      curr.p <- prop.p
      curr.vals <- prop.vals
      result[i,10] <- 1 	# indicate acceptance
    } else {
      result[i,10] <- 0 	# indicate rejection
    }
    
    result[i,1:6] <- c(curr.p,TE.imp)
  } # end of iter loop
  result
}


# simulate data from reflect Beta distribution
rrefbeta <- function(n,lambda)  {
  if(lambda<=0)  { 
    return(rbeta(n, 1, 1-lambda))
  }  else  return(rbeta(n, 1+lambda, 1))
}





















