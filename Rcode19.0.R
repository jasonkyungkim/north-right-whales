# ----------------------------------------------------------------------------- #

# ------- Estimating the populaton size of North Atlantic right whales ------- #

# R code by Steve Wang, Ling Zhong, Thibault Vernier, Claudia Xu, Tai Thongthai

# v1, 5/26/2011: first working version, but database not yet arranged into distinct sightings
# v2, 5/29/2011: adds code by Ling for processing database into distinct sightings
# v3, 5/31/2011: merges sightings that are close together in time
# v4, 5/31/2011: replaces rdist.earth with a faster function, gcd.slc
# v5, 6/02/2011: adds cutoff for distance traveled and check for sightingdb errors
# v6, 6/06/2011: adds shell call to CatchAll and samples M from negative binomial
# v7, 6/08/2011: saves CatchAll results in separate folders
# v8, 6/09/2011: samples from Poisson if SE ~ Mean; adds maxageclass
# v9, 6/16/2011: adds new version of merging rule, creates whale table
# v10, 6/23/2011: adds new adjustment for P(alive) based on behavior
# v11, 6/28/2011: fixed merging rule (old version wasnât keeping singletons)
# v12, 6/28/2011: updated whalelist to eliminate known dead whales
# v12b, 6/29/2011: finds coeffs for each sex separately (added by Ling)
# v12c, 7/4/2011: finds coeffs for each sex/age separately (added by Ling)
# v13, 7/05/2011: returns to finding coeffs for all whales pooled (not separately by sex/age)
#                 fixed gcd function to avoid NaN when whale doesnât move
#                 added code for plotting original vs updated P(alive)
# v14, 7/11/2011: removes sightings with date = 0
# v15, 8/20/2011: improves checks for whales going too fast
# v16, 8/05/2019: TV parallelized code and made other updates over summer 2019
# v17, 5/29/2020: SCW cleaned up code; put functions in separate file
# v18, 6/18/2020: JK and HS split covariates into 'ever' and 'now'
# v18.x, 9/25/2020: moved age cleaning before the regression
# v18.x, 12/2/2020: VDA added code that runs both Solow and ABM methods and compare results
# v19.0, 9.24.2021: SCW cleaned up code, changed variable names, changed indexing to subset(), etc.



# ----------------------------- Preliminaries -------------------------------- #

# clear workspace
rm(list=ls())

# load packages
library(foreach)
library(doParallel)
library(breakaway)
library(maps)
#library(CatchAll)
#library(tidyverse)
library(lubridate)
library(dplyr)
library(DHARMa)

# set working directory
setwd("C:/Users/bby0537/Documents/GitHub/north-right-whales")

# load functions
source("abm38.2inv.R")

# ---------------------------- Process behaviors ------------------------------ #
source("SemiAutomaticBehaviorParsing3.2.R")

# ------------------------ Read in sightings database ------------------------- #

# set initial and final year of sightingdb
source("whalefunctions.R")
initialyear <- 1986
finalyear <- 2017
now <- converttime(finalyear, 12, 31, 2359)

# read in sightings database; assumes database is sorted by whale, then date
sightingdb <- read.table("databasetemp.txt", header=T, stringsAsFactors=T) 

# remove sightings that are too new
sightingdb <- subset(sightingdb, sightingdb$SightingYear <= finalyear)
# don't delete sightings that are too old (before the initialyear) 
#  they are needed to calculate probalive for whales seen after initialyear
numsightings <- nrow(sightingdb)                      # number of sightings in database

# ---------------------- Read in vessel strike database ----------------------- #
strikedb <- read.csv('Vessel strike events for Steve Wang.csv')
# note: whales 1128, 1501, 1504 seen before 1986
# strikedb[order(strikedb$VesselStrikeId),]    # erase?
struckID <- strikedb[,"EGNo"]
numwhalesstruck <- length(unique(struckID))

# -------------------------- Create whale database ---------------------------- #

# create table of individual whales and their covariates
id <- sightingdb[,"SightingEGNo"]                     # whale ID number
idlist <- unique(id)                                  # list of distinct whale ID numbers
numwhales <- length(idlist)                           # number of sighted whales
whaledb <- matrix(NA,numwhales,29)                    # whale-level database
colnames(whaledb) <- c("id","sex","age","firstyear","lastyear","sightings",
                          "sagEver","sagNow","feedEver","feedNow","momwcalfNow","momwcalfEver", 
                          "entglEver","entglNow","fentglyear","disentglEver","disentglNow","disentglyear",
                          "yearsentgl","sickEver","sickNow","medicalEver","medicalNow","dead",
                          "strikes", "strike_severity", "deepStrike","lastGap","maxGap")

zero_days <- length(sightingdb["SightingDay"][sightingdb["SightingDay"]==0])
zero_months <- length(sightingdb["SightingMonth"][sightingdb["SightingMonth"]==0])

days_dist <- sightingdb$SightingDay[sightingdb$SightingDay != 0]
months_dist <- sightingdb$SightingMonth[sightingdb$SightingMonth != 0]

#replace 0 months and days with random values (from the actual distribution)
sightingdb["SightingMonth"][sightingdb["SightingMonth"]==0] <- sample(months_dist,size=zero_months,replace=TRUE)
sightingdb["SightingDay"][sightingdb["SightingDay"]==0] <- sample(days_dist,size=zero_days,replace=TRUE)

year <- sightingdb$SightingYear
month <- sightingdb$SightingMonth
day <- sightingdb$SightingDay

date <- as.Date(paste(year, "-", month, "-", day, sep=""))
sightingdb$sighting_full_date <- date

while (sum(is.na(sightingdb$sighting_full_date)) != 0){
  num_na<- length(which(is.na(sightingdb$sighting_full_date)))
  for (i in 1:num_na){
    year <- as.integer(sightingdb[which(is.na(sightingdb$sighting_full_date))[num_na],5])
    month <- as.integer(sightingdb[which(is.na(sightingdb$sighting_full_date))[num_na],6])
    day <- sample(days_dist,size=1,replace=TRUE)
    
    date <- as.Date(paste(year, "-", month, "-", day, sep=""))
    
    sightingdb[which(is.na(sightingdb$sighting_full_date))[num_na],25] <- date
  }
}




for(i in 1:numwhales)  {

  currid <- as.numeric(idlist[i])                     # id number of current whale
  subset <- (1:numsightings)[id==currid]              # currsightings corresponding to this whale
  currsightings <- sightingdb[subset,]                # submatrix of sightings of this whale
  k <- length(subset)                                 # number of sightings of this whale
  
  currsightings <- currsightings[order(currsightings$sighting_full_date),]
  
  if (k>1){
  gaps <- as.numeric((currsightings$sighting_full_date[2:k] - currsightings$sighting_full_date[1:k-1]))
  }
  
  # basic data
  whaledb[i,1] <- currid                              # id number of this whale
  whaledb[i,2] <- currsightings[,"GenderCode"][k]     # gender of this whale
  whaledb[i,3] <- currsightings[,"AgeClassCode"][k]   # age class as of latest sighting
  whaledb[i,4] <- min(currsightings[,"SightingYear"]) # first year sighted
  whaledb[i,5] <- max(currsightings[,"SightingYear"]) # last year sighted 
  whaledb[i,6] <- k                                   # how many sightings of this whale
  whaledb[i,7] <- max(currsightings[,"SAG"])          # ever seen in SAG
  whaledb[i,8] <- currsightings[,"SAG"][k]            # seen in SAG at latest sighting
  whaledb[i,9] <- max(currsightings[,"FEED"])         # ever sighted feeding
  whaledb[i,10] <- currsightings[,"FEED"][k]          # sighted feeding at latest sighting
  whaledb[i,11] <- currsightings[,"MOMWCALF"][k]      # mom has calf at latest sighting
  whaledb[i,12] <- max(currsightings[,"MOMHADCALF"])  # mom has ever had calf
  if(whaledb[i,2]==2)  whaledb[i,12] <- 0             # if a whale is male, set momhadcalf to 0
  whaledb[i,13] <- max(currsightings[,"EVERENTGL"])   # ever entangled
  whaledb[i,14] <- currsightings[,"NOWENTGL"][k]      # entangled at latest sighting
  whaledb[i,16] <- max(currsightings[,"DISENTGL"])    # ever disentangled
  if(whaledb[i,13]==0)   whaledb[i,16] <- 0           # DISENTGL = 0 if whale never entangled
  whaledb[i,17] <- currsightings[,"DISENTGL"][k]      # disentangled at latest sighting
  whaledb[i,20] <- max(currsightings[,"SICK"])        # ever sighted sick
  whaledb[i,21] <- currsightings[,"SICK"][k]          # latest sighting sick
  whaledb[i,22] <- max(currsightings[,"MEDICAL"])     # ever given medical assistance
  whaledb[i,23] <- currsightings[,"MEDICAL"][k]       # latest sighting given medical assistance
  whaledb[i,24] <- max(currsightings[,"DEAD"])        # known to be dead
  
  # entanglement data
  if(max(currsightings[,"FENTGL"])==1)  {             # if FENGTL, record the first year
    whaledb[i,15] <- min( subset(currsightings, FENTGL==1)$SightingYear )
  }  else  {
    whaledb[i,15] <- NA
  }
  if(max(currsightings[,"DISENTGL"])==1)  {           # if DISENGTL, record the last year
    whaledb[i,18] <- max( subset(currsightings, DISENTGL==1)$SightingYear )
  }  else  {
    whaledb[i,18] <- NA
  }  
  # calculate years entangled
  # if FENTGL but not DISENTGL, then assume still entangled
  # *** won't work if entangled >1 times ***
  if(!is.na(whaledb[i,15]) & is.na(whaledb[i,18]))  {   
    whaledb[i,19] <- finalyear - whaledb[i,15] + 1 
  }  else  {                 # has been disentangled
    whaledb[i,19] <- whaledb[i,18] - whaledb[i,15] + 1
  }

  # vessel strike data
  if(currid %in% struckID)  {					      # if this whale has been struck, get strike info
  	whaledb[i,25] <- max( subset(strikedb[,"EventNo"], struckID==currid) )      # number of strikes
    # save severity of strike
    if(max(grepl("SUPERFICIAL", subset(strikedb[,"VesselStrikeComment"], struckID==currid), fixed="True")))  {
      whaledb[i, 26] <- 1                             # severity of strike
      whaledb[i, 27] <- 0                             # deep strike?
    } 
    if(max(grepl("SHALLOW", strikedb[,"VesselStrikeComment"][struckID==currid], fixed="True")))  {
      whaledb[i, 26] <- 2                             # severity of strike
      whaledb[i, 27] <- 0                             # deep strike?
    } 
    if(max(grepl("DEEP", strikedb[,"VesselStrikeComment"][struckID==currid], fixed="True")))  {
      whaledb[i, 26] <- 3                             # severity of strike
      whaledb[i, 27] <- 1                             # deep strike?
    }
  }  else  {
  	whaledb[i,25] <- 0							      # number of strikes = 0 if whale has not been struck
  	whaledb[i,26] <- 0								  # strike_severity = 0 if whale has not been struck
  	whaledb[i,27] <- 0								  # deepStrike = 0 if whale has not been struck
  }
  
  gaps
  
  # last gap
  if(k > 1){
    whaledb[i,28] <- gaps[k-1]
    whaledb[i,29] <- max(gaps)
  }else{
    # last gap
    whaledb[i,28] <- gaps[1]
    # max gap
    whaledb[i,29] <- gaps[1]
  }

  # gaps <- NULL
  # note: if a whale was struck multiple times, the most severe strike type will be recorded
    
}

df_whaledb <- data.frame(whaledb)

# clean up and recode data
# note: we often use 0 instead of NA so that in the regression, inapplicable whales 
#   will have 0 adjustment rather than being NA and dropped from the regression

# remove whales whose sightings entirely predate initialyear
whaledb <- subset(whaledb, whaledb[,"lastyear"] >= initialyear)
numwhales <- nrow(whaledb)

# save whaledb
allwhaledb <- whaledb                        # save a copy in case needed later
write.table(whaledb, "whaledb.txt", sep="\t", row.names=F)



# ------------------------- Process sighting database -------------------------- #

# delete sightings whose whale ID numbers are not in the database (e.g., whales never sighted post-1986)
indices_to_delete <- c() 
for(i in 1: numsightings)  {
  if(!((sightingdb$SightingEGNo)[i] %in% whaledb[,"id"]))
    indices_to_delete <- c(indices_to_delete,i)
}
sightings_to_keep <- rep(TRUE, numsightings)
sightings_to_keep[indices_to_delete] <- FALSE
sightingdb <- subset(sightingdb, sightings_to_keep)
numsightings <- dim(sightingdb)[1]   # revised number of sightings, excluding whales not sighted post-1986


# assign maximum (latest) age class of each whale to all sightings of that whale
# this is used only for selecting subsets of the sightings database 
#   (e.g., selecting sightings of all current adults, even if they weren't adults at the time)
# alternately, could just take the last age recorded, but this way is safer if later ages were misrecorded  
tempid <- sightingdb[,"SightingEGNo"]                 # whale ID numbers currently; for convenience
maxageclass <- sightingdb[,"AgeClassCode"]            # initialize
for(i in (numsightings-1):1)  {                       # loop backwards over all sightings
  if(tempid[i]==tempid[i+1])                          # if this is the same whale
    maxageclass[i] <- maxageclass[i+1]                #   then use the later age class
  # otherwise, this is the latest sighting of a new whale; leave its maxageclass as is
}
sightingdb <- cbind(sightingdb, maxageclass)
rm(tempid, maxageclass)


# verify age class classifications
agenum <- sightingdb[,"Age"]
agenum = factor(agenum,levels(agenum)[c(1:3,14,25,33:38,4:13, 15:24, 26:32, 39,40)])
levels(agenum)
table(agenum, sightingdb[,"AgeClassCode"])



# ----------------- Create variables for logistic regression ------------------ #

# recode sex
female <- (whaledb[,"sex"])
female[female==2] <- -1                                    # -1 = males, 1 = females
# Changed from NA to 0 so ~50 whales would not be thrown out
female[female==3] <- 0 # NA                                # 0 = unknown

# recode age
# in whaledb, age is discretized as 1=adult, 2=calf, 3=juvenile, 4=unknown
# recode whales now at least 9 years old as an adult (even if they weren't when the data were collected)
for(i in 1:numwhales)  {
  minage <- finalyear - whaledb[,"firstyear"][i]
  if (minage >= 9)    whaledb[,"age"][i] <- 1              # 1 = adult
}
adult <- whaledb[,"age"]
adult[adult==2 | adult==3] <- -1                           # -1 = calf/juvenile
adult[adult==4] <- 0   # NA                                # 0 = unknown

# define and recode other variables
momwcalfNow <- whaledb[,"momwcalfNow"]                     # 1 = (adult female) w calf
momwcalfNow[momwcalfNow==0 & female==1 & adult==1] <- -1   # -1 = adult female w/o calf
                                                           # 0 = (male or juvenile)
momwcalfEver <- whaledb[,"momwcalfEver"]
momwcalfEver[momwcalfEver==0 & female==1 & adult==1] <- -1 # similar to above

firstyear <- whaledb[,"firstyear"]
lastyear <- whaledb[,"lastyear"]
sightings <- whaledb[,"sightings"]
sagEver <- whaledb[,"sagEver"]
sagNow <- whaledb[,"sagNow"]
feedEver <- whaledb[,"feedEver"]
feedNow <- whaledb[,"feedNow"]
entglEver <- whaledb[,"entglEver"]
entglNow <- whaledb[,"entglNow"]
fentglyear <- whaledb[,"fentglyear"]
disentglEver <- whaledb[,"disentglEver"]
disentglNow <- whaledb[,"disentglNow"]
disentglyear <- whaledb[,"disentglyear"]
yearsentgl <- whaledb[,"yearsentgl"]
sickNow <- whaledb[,"sickNow"]
sickEver <- whaledb[,"sickEver"]
sick <- sickNow | sickEver
medical <- whaledb[,"medicalEver"] | whaledb[,"medicalNow"]
sickmed <- sick | medical
dead <- whaledb[,"dead"]
alive <- 1 - dead    # alive means not confirmed dead; doesn't mean confirmed alive
strikes <- whaledb[,"strikes"]
strike_severity <- whaledb[,"strike_severity"]
strike_severity <- factor(strike_severity, ordered=TRUE)
strikeEver <- as.numeric(strikes > 0)
deepStrike <- whaledb[,"deepStrike"]

# recode disentglEver to separate out never-entangled (0) and not disentangled (-1)
disentglEver[entglEver==1 & disentglEver==0] <- -1

# add new variables to database; replace recoded variables
whaledb <- cbind(whaledb, adult, female, strikeEver)
whaledb[,"disentglEver"] <- disentglEver
whaledb[,"momwcalfEver"] <- momwcalfEver
whaledb[,"momwcalfNow"] <- momwcalfNow

# The above is potentially risky/confusing. Would be better to just refer to the 
#  variables as columns of whaledb, rather than defining new vectors, since we will
#  later edit whaledb, but these edits are not automatically propagated to the 
#  new vectors.


# ------------------------- Logistic regression models -------------------------- #

# all covarites except fentglyear, disentglyear, yearsentgl, 
# and excluding feedNow, sagNOW and disentglNow, which have huge SDs
# firstyear1 <- (firstyear-initialyear)
# firstyear2 <- (firstyear-initialyear)^2
# fit0 <- (glm(alive ~ female + entglEver + entglNow + momwcalfEver + momwcalfNow 
#              + feedEver  + sagEver + disentglEver  
#              + sickmed + adult + firstyear1 + firstyear2 + deepStrike, 
#              family="binomial"));   summary(fit0)   

# omit sickmed, not collected for whole time period
# fit1 <- (glm(alive ~ female + entglEver + entglNow + momwcalfEver + momwcalfNow 
#              + feedEver  + sagEver + disentglEver  
#              + adult + firstyear1 + firstyear2, 
#              family="binomial"));   summary(fit1)  

# omit sagEver and feedEver since most whales have been seen doing these behaviors
# fit2 <- (glm(alive ~ female + entglEver + entglNow + momwcalfEver + momwcalfNow 
#              + disentglEver  
#              + adult + firstyear1 + firstyear2, 
#              family="binomial"));   summary(fit2)  

# other possible models
# fit0 <- (glm(alive ~ female + entglEver + momwcalfEver + sagNow, family="binomial"));   summary(fit0)  
# fit1 <- (glm(alive ~ female + momwcalfEver + entglEver, family="binomial"));   summary(fit1)   
# fit2 <- (glm(alive ~ adult + momwcalfNow + momwcalfEver + female + entglEver, family="binomial"));   summary(fit2)  
# fit2 <- (glm(alive ~ adult + momwcalfEver + female + entglEver + entglNow, family="binomial"));   summary(fit2)  
# fit4 <- (glm(alive ~ adult + calf + female + entgl, family="binomial"));   summary(fit4)   
# fit4 <- (glm(alive ~ adult + momwcalfNow + momwcalfEver + feedEver + entglEver + entglNow, family="binomial"));   summary(fit4)   

# these models work fairly well
# fit3 <- (glm(alive ~ female + momwcalfEver + momwcalfNow + entglEver + entglNow + adult + disentglEver, family="binomial"));  summary(fit3) 
# fit3 <- (glm(alive ~ female + momwcalfEver + momwcalfNow + entglEver + entglNow + adult, family="binomial"));   summary(fit3)  
# fit3 <- (glm(alive ~ female + momwcalfEver + momwcalfNow + entglEver + entglNow + adult + disentglNow, family="binomial"));  summary(fit3)   
# disentglNow has huge SE
# female has small coefficient, could omit it

# add vessel strike info -- works much better. use this as final choice
LRfit <- (glm(alive ~ adult + female + feedEver + momwcalfEver + momwcalfNow + entglEver + entglNow 
            + disentglEver + strikeEver + deepStrike, family="binomial"));  summary(LRfit) 
coeffs <- coef(LRfit)[2:11]
coeffs

  
# check logistic regression model assumptions
# the list vector below should have all the models that are being considered
modelFits <- list(LRfit)
par(mfrow=c(2,2))
# # this figure represents a setting of the assumptions of the model
# for (i in 1:length(modelFits)) {
#   # simulate residuals in DHARMa
#   res.modelFits <- simulateResiduals(modelFits[[i]], refit = F, n = 1000)
#   
#   # check for uniformity
#   testUniformity(res.modelFits)
#   
#   # test for overdispersion
#   testDispersion(res.modelFits) #this uses the dharma simulaiton output, notice the res. before models
#   
#   # test for zero inflation
#   testZeroInflation(res.modelFits) #this uses the dharma simulaiton output, notice the res. before models
#   
#   # plot the resiuals of each variable - so we are going to plot Population
#   plotResiduals(res.modelFits, asFactor = T)
#   # plotResiduals(res.modelFits, asFactor = F)
#   
#   # VDA: none of the populations are too terrible in residuals so we could keep them
#   summary(modelFits[[i]])
# }


# ---------------------------- Further processing ----------------------------- #

# note this cleaning code must come after the whaledb is created, because we don't want to 
#   remove sightings that might contain information about covariates before the logistic regression



# ---------------------------- Remove dead whales ----------------------------- #

# remove confirmed dead whales from sighting and whale databases
# this processing must occur after logistic regression, since L.R. must include dead whales
# note: whale 1128 was already removed since all its sightings predate 1986
deadID <- subset(sightingdb[,"SightingEGNo"], sightingdb[,"DEAD"]==1)    # IDs of dead whales
deadID <- unique(deadID)                        # remove duplicate ID numbers
for(i in 1:length(deadID))                      # remove dead whales from sighting database
  sightingdb <- subset(sightingdb, sightingdb[,"SightingEGNo"] != deadID[i])
whaledb <- subset(whaledb, dead==0)             # remove dead whales from whale database

# don't delete; 0 day to random number between 1 and 31 and 0 month from actual distribution of months; sort; then calculate
# remove sightings with month = 0 or day = 0  (only a few dozen of these; just delete them)
sightingdb <- subset(sightingdb, sightingdb[,"SightingDay"] != 0)
# should have been removed earlier (< 1986)
sightingdb <- subset(sightingdb, sightingdb[,"SightingMonth"] != 0)   

# update database info
id <- sightingdb[,"SightingEGNo"]               # whale ID number
idlist <- unique(id)                            # list of distinct whale ID numbers
numwhales <- length(idlist)                     # number of sighted whales
numsightings <- nrow(sightingdb)                # number of sightings in database

# ------ Remove whales travelling unrealistically fast between sightings ------- #
# ------------------ and multiple non-independent sightings -------------------- #

# These are done together because they both depend on gaps between sightings (waittime).
# ***** Not clear that too-fast sightings should be deleted. Even if location is wrong, they are actual sightings.

# get sighting date and time
year <- sightingdb[,"SightingYear"]                        # year of sighting
month <- sightingdb[,"SightingMonth"]                      # month of sighting
day <- sightingdb[,"SightingDay"]                          # day of sighting
time <- sightingdb[,"SightingTime"]                        # time of sighting

# impute values for missing times (where time is recorded as 0 (i.e., midnight))
nonzerotimes <- subset(time, time != 0)
nummissingtimes <- length(which(time==0))                  # 1742 sightings have missing times
# qqnorm(nonzerotimes)                                     # dist. of all times looks reasonably normal
# assign all missing times to be noon
imputedtimes <- rep(1200, nummissingtimes)
# not used -- impute random times
#imputedtimes <- rnorm(nummissingtimes, mean(nonzerotimes), sd(nonzerotimes))
#imputedtimes[imputedtimes < 0] <- 1                       # check for impossible values
#imputedtimes[imputedtimes > 2359] <- 2359
#for(i in 1:nummissingtimes)                               # fix times that end in >59 minutes
#  if (imputedtimes[i] %% 100 > 59)    imputedtimes[i] <- imputedtimes[i] - 60 + 100
timeimp <- time                                            # initialize
timeimp[timeimp==0] <- round(imputedtimes)                 # replace 0 times with imputed times

# convert sighting dates and times to standard format (days since 1/1/0000)
stdtime <- converttime(year,month,day,timeimp) # convert to standard time format
sightingdb <- cbind(sightingdb, stdtime)                   # add standardized sighting time to database
rm(year, month, day, time, timeimp, imputedtimes, nonzerotimes)

# calculate waiting times in days (sighting gaps) since last sighting
waittime <- stdtime[2:numsightings] - stdtime[1:(numsightings-1)]   # days between consecutive sightings
waittime <- c(NA, waittime)                # first sighting in database has no waiting time
for(i in 2:numsightings)                   # if not the same whale, then make its waittime = NA
  if (id[i] != id[i-1])    waittime[i] <- NA

# check that we have the right number of NA values 
sum(is.na(waittime)) == numwhales          # should be TRUE
# imputing times causes a small number of sightings to be out of order (76 out of 63016)
#  arbitrarily set these waiting times to zero
waittime[waittime < 0] <- 0
# check to make sure there are no negative values (confirms that database is sorted correctly)
summary(waittime)
sightingdb$waittime <- round(waittime,5)   # add waiting time to sighting database

# calculate distance traveled (in km) since last sighting
lat <- sightingdb[,"Latitude"]                             # latitude of sighting
long <- sightingdb[,"Longitude"]                           # longitude of sighting
lat[lat==0] <- NA;    long[long==0] <- NA;                 # missing values were recorded as 0; change to NA
traveldist <- rep(NA, numsightings)
for(i in 2:numsightings)  {
  if (id[i] == id[i-1])  {                                 # if this is the same whale as previous row
    lat1 <- deg2rad(lat[i-1]);    long1 <- deg2rad(long[i-1])
    lat2 <- deg2rad(lat[i]);      long2 <- deg2rad(long[i])
    traveldist[i] <- gcd.slc(long1, lat1, long2, lat2)
  }
  # if not the same whale, leave its traveldist as NA
}
sightingdb$traveldist <- round(traveldist,5)               # add traveled distance to sighting database

# delete sightings whose traveldist/waittime combination exceeds top speed of whale
# delete multiple sightings within 30 days in the same region
maxspeed <- 9.3                                            # km/hr, from wikipedia
waittimecutoff <- 30                                       # days
sightingdb$speed <- round(traveldist/(waittime*24),3)
toofast <- rep(0, numsightings)
rowstokeep <- rep(1, numsightings)
region <- sightingdb[,"RegionCode"]

for(i in 2:numsightings)  {
  if(!is.na(traveldist[i]) & !is.na(waittime[i]) & traveldist[i] > waittime[i]*9.3*24)  {
    rowstokeep[i] <- 0                                     # exceeds top speed of whale
    toofast[i] <- 1
  }
  if( !is.na(waittime[i]) & waittime[i]<=waittimecutoff & region[i]==region[i-1] )
    rowstokeep[i] <- 0                                     # multiple sightings in same area within short time
}

# save list of too-fast whales before deleting them from the database
write.table(sightingdb[toofast==1,], "toofast.txt", sep="\t", row.names=F)
sightingdb <- sightingdb[rowstokeep==1,]

# update database info
id <- sightingdb[,"SightingEGNo"]               # whale ID number
idlist <- unique(id)                            # list of distinct whale ID numbers
numwhales <- length(idlist)                     # number of sighted whales
numsightings <- nrow(sightingdb)                # number of sightings in database
rm(stdtime, traveldist, waittime)               # clean up

# --------------------------- Optional: select subsets ------------------------ #

# select by sex
sex <- sightingdb[,"GenderCode"]
# sightingdb <- sightingdb[sex=="F",]

# select by age
maxageclass <- sightingdb[,"maxageclass"]
# sightingdb <- sightingdb[maxageclass=="A",]

# ------------- Calculate probabilities that each whale is alive -------------- #

# Set up parallel computing
cores <- detectCores()
cl <- makeCluster(cores) 
registerDoParallel(cl)
stdnow <- converttime(finalyear,12,31,2359)

# sdb <- data.frame(whaledb)

# Calculate P(alive) using Solow method
p_solow <- foreach(i = 1:numwhales, .combine="c", .inorder=TRUE)  %dopar%  {
  currid <- whaledb[i,"id"]                            # ID number of current whale
  currsightings <- sightingdb$stdtime[id==currid]      # select matching sightings
  
  if(df_whaledb$sightings[df_whaledb$id==currid] < 4){
    k <- df_whaledb$lastGap[df_whaledb$id==currid]
    
    a <- nrow(subset(df_whaledb, lastGap > k & dead == 0))
    
    b <- nrow(subset(df_whaledb, lastGap > k | maxGap > k))
    
    temp <- a/b
    
  }else{
    temp <- solow(currsightings, now=stdnow)
  }
  temp
}
# Calculate P(alive) using ABM method -- giving errors due to numerical instability, skip for now
#p_abm <- foreach(i = 1:numwhales, .combine="c", .inorder=TRUE)  %dopar%  {
#  currid <- whaledb[i,"id"]                            # ID number of current whale
#  # for ABM, use ext=T, dist=T
#  currsightings <- sightingdb$stdtime[id==currid]      # select matching sightings
#  temp <- 1 - abm38inv(currsightings, ext=T, distance=T, now=stdnow)
#  temp
#}
stopCluster(cl)
p_abm <- rep(1, numwhales)    # temp values to make sure code doesn't crash

# Add to whale database
p_solow <- round(p_solow, 4)
p_abm <- round(p_abm, 4)
whaledb <- cbind(whaledb, p_solow, p_abm)

# Plot results
pdf("whaleplots.pdf")
par(mfrow=c(2,2))
hist(p_solow, breaks=30, main=paste("sum =", round(sum(p_solow))))
hist(p_abm, breaks=30, main=paste("sum =", round(sum(p_abm))))
plot(p_solow, p_abm)
plot(p_abm, p_solow)
boxplot(p_solow ~ whaledb[,5], main="p_solow by last sighting year")
boxplot(p_abm ~ whaledb[,5], main="p_abm by last sighting year")
# draw map of sightings
lat <- sightingdb[,"Latitude"];    long <- sightingdb[,"Longitude"]
plot(lat ~ long, pch=1, cex=.5, col="red")
maps::map(add=T)

# ------------- Adjust probabilities based on covariates -------------- #

# Convert probabilities to odds
odds_solow <- p_solow/(1-p_solow)
odds_abm <- p_abm/(1-p_abm)

# Adjust solow probabilities of being alive
coeffs
for(i in 1:numwhales)  {
  currid <- whaledb[i,"id"]                             # ID number of current whale
  currwhale <- whaledb[whaledb[,"id"]==currid,]         # covariates for this whale
  # choose the covariates we used in the logistic regression above (must be in correct order)
  x <- currwhale[c("adult", "female", "feedEver", "momwcalfEver", "momwcalfNow", 
                   "entglEver", "entglNow", "disentglEver", "strikeEver", "deepStrike")] 
  adjustment <-  exp(sum(coeffs*x))
  odds_solow[i] <- odds_solow[i] * adjustment
}
p_solow_adj <- odds_solow/(1 + odds_solow)
print(sum(p_solow_adj))

# Adjust abm probability of being alive
for(i in 1:numwhales)  {
  currid <- whaledb[i,"id"]                             # ID number of current whale
  currwhale <- whaledb[whaledb[,"id"]==currid,]         # covariates for this whale
  # choose the covariates we used in the logistic regression above (must be in correct order)
  x <- currwhale[c("adult", "female", "feedEver", "momwcalfEver", "momwcalfNow", 
                   "entglEver", "entglNow", "disentglEver", "strikeEver", "deepStrike")] 
  adjustment <-  exp(sum(coeffs*x))
  odds_abm[i] <- odds_abm[i] * adjustment
}
p_abm_adj <- odds_abm/(1 + odds_abm)
# print(sum(p_abm_adj))

# Check for whales with p_alive = 1, which will give Inf odds
p_solow_adj[p_solow==1] <- 1

# whale 3030 returns p_abm = 1 and therefore odds = Inf and adjusted value = NaN.
#   This is a numerical error, as this whale is likely dead.

# Plot old and new probability of survival
par(mfrow=c(2,2))
hist(p_solow, breaks = 30, xlim = c(0,1), main=paste("sum =", round(sum(p_solow))))
hist(p_solow_adj, breaks = 30, xlim = c(0,1), main=paste("sum =", round(sum(p_solow_adj))))
#hist(p_abm, breaks = 30, xlim = c(0,1), main=paste("sum =", round(sum(p_abm))))
#hist(p_abm_adj, breaks = 30, xlim = c(0,1), main=paste("sum =", round(sum(p_abm_adj))))

# Make plots of new vs original p_alive values by covariates (solow only)

plot(p_solow_adj ~ p_solow, type="n", xlim = c(0,1), ylim = c(0,1))
points(p_solow_adj ~ p_solow, col="black", cex=.7)
abline(0,1, col="gray")
title(main="all")

plot(p_solow_adj ~ p_solow, type="n", xlim = c(0,1), ylim = c(0,1))
subset <- whaledb[,"female"] == 1
 points(p_solow_adj[subset] ~ p_solow[subset], col="red", cex=.7)
subset <- whaledb[,"female"] == -1
 points(p_solow_adj[subset] ~ p_solow[subset], col="blue2", cex=.7)
subset <- whaledb[,"female"] == 0
 points(p_solow_adj[subset] ~ p_solow[subset], col="gray", cex=.7)
abline(0,1, col="gray")
title(main="sex")
legend(-.05,1.05, c("female","male","unknown"), fill=c("red", "blue2","gray"), cex=.6, border=F, bty="n")

plot(p_solow_adj ~ p_solow, type="n", xlim = c(0,1), ylim = c(0,1))
subset <- whaledb[,"strikeEver"] == 1
 points(p_solow_adj[subset] ~ p_solow[subset], col="red", cex=.7)
subset <- whaledb[,"strikeEver"] == 0
 points(p_solow_adj[subset] ~ p_solow[subset], col="blue2", cex=.7)
abline(0,1, col="gray")
title(main="strike ever")
legend(-.05,1.05, c("struck","not struck"), fill=c("red", "blue2"), cex=.6, border=F, bty="n")

plot(p_solow_adj ~ p_solow, type="n", xlim = c(0,1), ylim = c(0,1))
subset <- whaledb[,"deepStrike"] == 1
 points(p_solow_adj[subset] ~ p_solow[subset], col="red", cex=.7)
subset <- whaledb[,"deepStrike"] == 0
 points(p_solow_adj[subset] ~ p_solow[subset], col="blue2", cex=.7)
abline(0,1, col="gray")
title(main="deep strike")
legend(-.05,1.05, c("deep strike","not deep strike"), fill=c("red", "blue2"), cex=.6, border=F, bty="n")

plot(p_solow_adj ~ p_solow, type="n", xlim = c(0,1), ylim = c(0,1))
subset <- whaledb[,"momwcalfNow"] == 1
 points(p_solow_adj[subset] ~ p_solow[subset], col="red", cex=.7)
subset <- whaledb[,"momwcalfNow"] == 0
 points(p_solow_adj[subset] ~ p_solow[subset], col="darkgray", cex=.7)
subset <- whaledb[,"momwcalfNow"] == -1
 points(p_solow_adj[subset] ~ p_solow[subset], col="blue2", cex=.7)
abline(0,1, col="gray")
title(main="with calf")
legend(-.05,1.05, c("w/ calf","NA","w/o calf"), fill=c("red", "darkgray","blue2"), cex=.6, border=F, bty="n")

plot(p_solow_adj ~ p_solow, type="n", xlim = c(0,1), ylim = c(0,1))
subset <- whaledb[,"entglEver"]==1 & whaledb[,"entglNow"]==1
 points(p_solow_adj[subset] ~ p_solow[subset], col="red", cex=.7)
subset <- whaledb[,"entglEver"]==0
 points(p_solow_adj[subset] ~ p_solow[subset], col="blue2", cex=.7)
subset <- whaledb[,"entglEver"]==1 & whaledb[,"entglNow"]==0
 points(p_solow_adj[subset] ~ p_solow[subset], col="purple", cex=.7)
abline(0,1, col="gray")
legend(-.05,1.05, c("entangled now","never entangled","ever entangled, not now"), fill=c("red","blue2", "purple"), cex=.6, border=F, bty="n")
title(main="entangled")

plot(p_solow_adj ~ p_solow, type="n", xlim = c(0,1), ylim = c(0,1))
subset <- whaledb[,"feedEver"]==1
 points(p_solow_adj[subset] ~ p_solow[subset], col="red", cex=.7)
subset <- whaledb[,"feedEver"]==0
 points(p_solow_adj[subset] ~ p_solow[subset], col="blue2", cex=.7)
abline(0,1, col="gray")
legend(-.05,1.05, c("seen feeding","not seen feeding"), fill=c("red","blue2"), cex=.6, border=F, bty="n")
title(main="feeding")

plot(p_solow_adj ~ p_solow, type="n", xlim = c(0,1), ylim = c(0,1))
subset <- subset <- whaledb[,"age"]==1
 points(p_solow_adj[subset] ~ p_solow[subset], col="blue2", cex=.7)
subset <- subset <- whaledb[,"age"]==2  |   whaledb[,"age"]==3
 points(p_solow_adj[subset] ~ p_solow[subset], col="red", cex=.7)
subset <- subset <- whaledb[,"age"]==4
 points(p_solow_adj[subset] ~ p_solow[subset], col="gray", cex=.7)
abline(0,1, col="gray")
legend(-.05,1.05, c("adult","calf/juvenile", "unknown"), fill=c("blue2","red","gray"), cex=.7, border=F, bty="n")
title(main="age")



# ------------------------- Estimate number of unseen whales --------------------------- #

# Initialize
numreps <- 10                                       # number of random samples to draw
cores <- detectCores()                              # set up parallel computation
cl <- makeCluster(cores) 
registerDoParallel(cl)

# Main loop
results_solow <- foreach(rep=1:numreps, .combine=rbind, .inorder=FALSE) %dopar% {
  
  cat(rep, " ")
  currsightings <- sightingdb                       # whales simulated as alive in this iteration
  theta <- rbinom(numwhales, 1, p_solow_adj)        # vector of 0 (dead) or 1 (alive) for each whale
  S <- sum(theta)                                   # number of living whales in this iteration
  
  # select only whales simulated as alive in this iteration
  currdead <- idlist[ (1:numwhales)[theta==0] ]     # get list of whales simulated as dead in this iteration 
  for(i in 1:length(currdead))  
    currsightings  <- currsightings[currsightings[,"SightingEGNo"] != currdead[i], ]   # keep omitting "dead" whales
  currid <- currsightings[,"SightingEGNo"]          # update ID numbers of remaining whales
  
  # create frequency of frequencies for use by breakaway
  freqs <- table(currid)                            # get frequency of sightings for each whale
  currtable <- table(freqs)                         # get frequency of frequencies
  input <- cbind(as.numeric(rownames(currtable)), currtable)    # convert to matrix
  colnames(input) <- c("", "f")

  # what is this?
  # input <- sightingdb.frame(index=input[,1],frequency=input[,2])    # Fixing sightingdb type for CatchAll input

  # invoke breakaway
  breakawayresults <- breakaway::breakaway(input)
  Nhat <- breakawayresults$est
  SE <- breakawayresults$seest
  model <- breakawayresults$name
  tau <- NA         # not sure what this is
  Mhat <- Nhat - S
  phihat <- NA
  alphahat <- NA
  
  if(!is.na(Nhat))  {     # sample from the appropriate distribution for the number of unsighted whales
    if(SE^2 <= Mhat)  {                           # sample from Poisson
      M <- rpois(1, Mhat)
    }
    if(SE^2 > Mhat)  {                            # overdispersed - sample from Negative Binomial
      # method of moments:
      # set the parameters of the negative binomial to the observed equivalents found by breakaway
      phihat <- Mhat/(SE^2)
      alphahat <- Mhat*phihat/(1-phihat)
      # draw the number of unsighted whales from the negative binomial
      M <- rnbinom(1, alphahat, phihat)
    }
    N <- M + S                                     # total number of whales = unsighted + sighted
  }
  
  c(tau, S, round(Nhat,3), round(SE,3), round(Mhat,3), round(alphahat,3), round(phihat,3), M, N, model)
}
stopCluster(cl)
colnames(results_solow) <- c("tau", "S", "Nhat", "SE", "Mhat", "alphahat", "phihat", "M", "N", "model")
head(results_solow)

# Get Nhat results as numeric
Nhat_results_solow <- as.numeric(results_solow[,"N"])

# make histogram of predicted N-hat values
write.table(results_solow, "results_solow.txt", sep="\t", row.names=FALSE)

par(mfrow=c(1,1))
hist(Nhat_results_solow, col="darkgray", border="white", ylab="", xlab="population size",
     main="Estimated population size", breaks = 20)

cat("mean of catchall predictions Nhat")
mean(na.omit(Nhat_results_solow))
abline(v=mean(na.omit(Nhat_results_solow)), col="red")
 
dev.off()

save.image("whaleresults.Rdata")

