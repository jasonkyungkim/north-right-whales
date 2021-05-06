# ---------------------------------------------------------------------------- #

# R code for Adaptive Beta Method
# version 38.2, July 30, 2019
# from the paper
# Adaptive credible intervals on stratigraphic ranges 
#   when recovery potential is unknown
# by Steve Wang, Phil Everson, Heather Zhou, David Chudzicki, Dasol Park
# Paleobiology 2015
# Before using, please contact the authors for updates, bug fixes, etc.:
# Steve Wang, scwang@swarthmore.edu

# ---------------------------------------------------------------------------- #




# --------------------------- Adaptive Beta method --------------------------- #

abm38inv <- function(x, distance=TRUE, ext=TRUE, base=NULL, now, PLOT=0)    {
  
  # returns the probability that the taxon is extinct at time/distance 'now'
  # x = vector of locations or dates of fossil occurrences
  # distance = T if measurements represent distance above a base,
  #            F if measurements represent time
  # ext = T if extinction, F if origination 
  # base = value to consider 0 if known; otherwise the min is used (if upperbound) and
  #         sample size is decreased by 1
  # now = the specific year we are interested in
  # PLOT = 1 to show plots
  

  # ---------- CHECK SAMPLE SIZE ---------- #
  n <- length(x)
  # for small samples, make up some arbitrary values
  if(n==1) return(.2)
  else if(n==2) return(.3)
  else if(n==3) return(.4)
  else if(n==4) return(.5)

  # if n ≥ 5, then use the ABM method
  else if(n>=5)  {
  
  
  # ---------- FUNCTION DEFINITIONS ---------- #
  
  # following are 'wrapper' functions required for R's integrate() function to work
  # note: these do not check for xmax < th; if xmax > th, an error will not be generated
  
  # prior mean and SD
  prmean <- 0
  prSD <- 2
  
  integrand.neglambdas <- function(L,th,x)  {
    # L = vector of lambdas,  th = theta value (scalar),  x = vector of strat. positions
    k <- length(L)    
    output <- rep(NA, k)
    for(i in 1:k) 
      output[i] <- ( sum( log( (1-L[i])/th * 1/(1-x/th)^L[i] ) ) )
    output <- output + dnorm(L, prmean,prSD, log=TRUE) + log(1/th)
    output <- exp(output)
    return(output)
  }
  
  integrand.poslambdas <- function(L,th,x)  {
    # L = vector of lambdas,  th = theta value (scalar),  x = vector of strat. positions
    k <- length(L)    
    output <- rep(NA, k)
    for(i in 1:k) 
      output[i] <- ( sum( log( (1+L[i])/th * (x/th)^L[i] ) ) )
    output <- output + dnorm(L, prmean,prSD, log=TRUE) + log(1/th)
    output <- exp(output)
    return(output)
  }
  
  integrand.thetasnegL <- function(th,L,x)  {
    # th = vector of thetas,  L = lambda value (scalar), x = vector of strat. positions
    k <- length(th) 
    output <- rep(NA, k)
    for(i in 1:k) 
      output[i] <- ( sum( log( (1-L)/th[i] * 1/(1-x/th[i])^L ) ) )
    output <- output + dnorm(L, prmean,prSD, log=TRUE) + log(1/th)
    output <- exp(output)
    return(output)  
  }
  
  integrand.thetasposL <- function(th,L,x)  {
    # th = vector of thetas,  L = lambda value (scalar), x = vector of strat. positions
    k <- length(th) 
    output <- rep(NA, k)
    for(i in 1:k) 
      output[i] <- ( sum( log( (1+L)/th[i] * (x/th[i])^L) ) )
    output <- output + dnorm(L, prmean,prSD, log=TRUE) + log(1/th)
    output <- exp(output)
    return(output)  
  }
  
  
  # pdf of reflected beta density (used in likelihood calculations)
  drefbeta <- function(x,L)  {
    if(L<=0)  { 
      return(dbeta(x, 1, 1-L))
    }  else  return(dbeta(x, 1+L, 1))
  }
  
  
  
  # ---------- PRE-PROCESS DATA ---------- #
  
  # If a base is specified, check that it is valid
  if(!is.null(base)) 
    if( (distance & ext)   & (base > min(x))  | 
        (distance & !ext)  & (base < max(x))  |
        (!distance & ext)  & (base < max(x))  |
        (!distance & !ext) & (base > min(x)) )
      stop("Invalid value for base")
  
  # Convert units relative to base of section or other zero point
  xraw <- x                            # save a copy for later use
  if((distance & ext) | (!distance & !ext))  {
    if(!is.null(base))                 # if a base is specified
      x <- x - base 
    if(is.null(base))  {               # if a base is not specified
      base <- min(x)                   # condition on smallest value and
      x <- x - base 
      x <- sort(x, decreasing=F)[-1]   #   reduce sample size by 1
    }
    
    # Scale 'now'
    now <- now - base
  } 
  
  if((distance & !ext) | (!distance & ext))  {
    if(!is.null(base))                 # if a base is specified
      x <- base - x
    if(is.null(base))  {               # if a base is not specified
      base <- max(x)                   # condition on smallest value and
      x <- base - x   
      x <- sort(x, decreasing=F)[-1]   #   reduce sample size by 1
    } 
  } 
  
  # Scale data so theta is approx. 100 (solely for numerical stability)
  xmax <- max(x)     
  simplethhat <- (n+1)/n * xmax
  scalefactor <- 100/simplethhat
  x <- x * scalefactor
  xmax <- max(x)
  
  # Scale 'now'
  now <- now*scalefactor
    
  # Set iteration parameters
  upperlimth <- 700;     numstepsth <- 1000
  lowerlimL <- -10;      upperlimL <-  10;      numstepsL <- 40
  Lvals <- seq(lowerlimL, upperlimL, length.out=numstepsL)
  Ldens <- rep(NA, numstepsL) 
  thetavals <- seq(xmax, upperlimth, length.out=numstepsth)
  thdens <- rep(NA, numstepsth)
  
  
  # ---------- ESTIMATE LAMBDA ---------- #
  
  # increment lambda values, integrating over theta values for each
  for(i in 1:numstepsL)  {
    Ldens[i] <- ifelse( Lvals[i]<=0, 
                        integrate(integrand.thetasnegL, xmax,upperlimth, L=Lvals[i], x=x)$value,
                        integrate(integrand.thetasposL, xmax,upperlimth, L=Lvals[i], x=x)$value )
  }
  # normalize lambda pdf to unit area  
  Ldens <- Ldens/sum(Ldens) 
  
  # calculate posterior quantities
  # cutoff <- which.max( cumsum(Ldens) >= .5 )       
  # Lmed <- Lvals[cutoff]                            # posterior median
  Lmean <- sum(Lvals*Ldens)                          # posterior mean
  Lhat <- Lmean                                      # point estimate for lambda
  Lvar <- sum(Ldens * (Lvals - Lmean)^2)             # posterior variance 
  # print(Lhat)                                      # for testing
  
  
  # ---------- ESTIMATE THETA ---------- #
  
  # increment theta values, integrating over lambda values for each
  for(i in 1:numstepsth)  {
    thdens[i] <- ( integrate(integrand.neglambdas, -Inf,0, th=thetavals[i], x=x)$value
                   + integrate(integrand.poslambdas,  0,Inf, th=thetavals[i], x=x)$value )
  }
  # normalize theta pdf to unit area  
  thdens <- thdens/sum(thdens) 
  

  # ---------- RETURN P(EXTINCT) ---------- #
  
  # find the index which 'now' corresponds to
  index <- which.max(thetavals >= now)

  # sum up the area under curve from base to now
  prob <- cumsum(thdens)
  
  return(prob[index])

}  # else n>=5
}



# testing data
# Anabarella data from Wang et al 2016 paper
x <- c(533.06579, 531.07032, 530.05212, 530.02947, 529.78322, 529.47184, 
       528.11654, 527.84072, 527.72420, 527.68698, 527.62881, 525.62909, 
       525.0029, 524.67879, 523.8070, 523.76623, 523.67817, 522.95226, 522.19971)

# Fake test data
y <- c(2010, 2012, 2012, 2014, 2015, 2017, 2017, 2019)

abm38inv(x, now=535.1)   # should be .5
abm38inv(x, now=542.4)   # should be .9


