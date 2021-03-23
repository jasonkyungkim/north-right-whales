# ----------------------------------------------------------------------------- #

# ------- Estimating the populaton size of North Atlantic right whales ------- #

# R code by Steve Wang, Ling Zhong, Thibault Vernier, Claudia Xu, Tai Thongthai

# v1, 5/26/2011: first working version, but database not yet arranged into distinct sightings
# v2, 5/29/2011: adds code by Ling for processing database into distinct sightings
# v3, 5/31/2011: merges sightings that are close together in time
# v4, 5/31/2011: replaces rdist.earth with a faster function, gcd.slc
# v5, 6/02/2011: adds cutoff for distance traveled and check for data errors
# v6, 6/06/2011: adds shell call to CatchAll and samples M from negative binomial
# v7, 6/08/2011: saves CatchAll results in separate folders
# v8, 6/09/2011: samples from Poisson if SE ~ Mean; adds maxageclass
# v9, 6/16/2011: adds new version of merging rule, creates whale table
# v10, 6/23/2011: adds new adjustment for P(alive) based on behavior
# v11, 6/28/2011: fixed merging rule (old version wasn’t keeping singletons)
# v12, 6/28/2011: updated whalelist to eliminate known dead whales
# v12b, 6/29/2011: finds coeffs for each sex separately (added by Ling)
# v12c, 7/4/2011: finds coeffs for each sex/age separately (added by Ling)
# v13, 7/05/2011: returns to finding coeffs for all whales pooled (not separately by sex/age)
#                 fixed gcd function to avoid NaN when whale doesn’t move
#                 added code for plotting original vs updated P(alive)
# v14, 7/11/2011: removes sightings with date = 0
# v15, 8/20/2011: improves checks for whales going too fast
# v16, 8/05/2019: TV parallelized code and made other updates over summer 2019
# v17, 5/29/2020: SCW cleaned up code; put functions in separate file
# v18, 6/18/2020: JK and HS split covariates into 'ever' and 'now'
# v19, 9/25/2020: moved age cleaning before the regression
# v20, 12/2/2020: VDA: added code that runs both Solow and ABM methods. Also added code
#                 to compare the difference between both methods because they show
#                 vastly different results

# ----------------------------- Preliminaries -------------------------------- #

# clear workspace
rm(list=ls())
dev.off()

# load packages
library(foreach)
library(doParallel)
library(breakaway)
#library(CatchAll)
library(tidyverse)
library(DHARMa)

# set working directory
# setwd("~/Documents/Research/Right whales/")
setwd("~/Desktop/2020 Summer Research")

# load database
source("SemiAutomaticBehaviorParsing3.2.R")

# load functions
#setwd("~/Documents/Research/Right whales/")
#setwd("/Users/vitor/Desktop/Right Whales/")
source("whalefunctions.R")
source("abm38.2inv.R")



# --------------------------- Pre-process dataset ---------------------------- #


# set working directory and read in entire database
# assumes data are sorted by whale, then date
data <- read.table("databasetemp.txt", header=T, stringsAsFactors=T)     # sightings-level database

# set initial and final year of data
initialyear <- 1986
finalyear <- 2017

# remove sightings that are too old or too new
# data <- data[data[,"SightingYear"]>=initialyear,]  # keeps only sightings post-1986; don't do this
data <- data[data[,"SightingYear"]<=finalyear,]
now <- as.numeric(converttime(finalyear, 12, 31, 2359))

#read in vessel strike data
vesselData <- read.csv('Vessel strike events for Steve Wang.csv')
#1128, 1501, 1504 seen before 1986

#head(vesselData)

#vesselData[order(vesselData$VesselStrikeId),]

vesselID <- vesselData[,"EGNo"]
nwhalestriked <- length(table(vesselID))

# make table of individual whales and their covariates
id <- data[,"SightingEGNo"]                              # whale ID number
nwhales <- length(table(id))                             # number of sighted whales
idlist <- rownames(table(id))                            # list of distinct whale ID numbers
numrows <- dim(data)[1]                                  # number of sightings in database
whaletable <- matrix(NA,nwhales,27)                      # whale-level database

for(i in 1:nwhales)  {
  tempid <- as.numeric(idlist[i])                        # id numberof current whale
  subset <- (1:numrows)[id==tempid]                      # rows corresponding to this whale
  sightings <- data[subset,]                             # submatrix of sightings of this whale
  k <- length(subset)                                    # number of sightings of this whale
  whaletable[i,1] <- tempid                              # id number of this whale
  whaletable[i,2] <- sightings[,"GenderCode"][k]         # gender of this whale
  whaletable[i,3] <- sightings[,"AgeClassCode"][k]       # age class as of latest sighting
  whaletable[i,4] <- sightings[,"SightingYear"][1]       # first year sighted
  whaletable[i,5] <- max(sightings[,"SightingYear"])     # last sighting year
  whaletable[i,6] <- k                                   # how many sightings of this whale
  whaletable[i,7] <- max(sightings[,"SAG"])              # ever seen in SAG
  whaletable[i,8] <- sightings[,"SAG"][k]                # latest SAG sighting
  whaletable[i,9] <- max(sightings[,"FEED"])             # ever sighted feeding
  whaletable[i,10] <- sightings[,"FEED"][k]              # sighted feeding at latest sighting
  whaletable[i,11] <- sightings[,"MOMWCALF"][k]          # mom has calf at latest sighting
  whaletable[i,12] <- max(sightings[,"MOMHADCALF"])      # mom has ever had calf
  whaletable[i,13] <- max(sightings[,"EVERENTGL"])       # ever entangled
  whaletable[i,14] <- sightings[,"NOWENTGL"][k]          # latest entanglement
  whaletable[i,16] <- max(sightings[,"DISENTGL"])        # ever disentangled
  if(whaletable[i,13]==0)   whaletable[i,16] <- 0        # DISENTGL = 0 if whale never entangled
  whaletable[i,17] <- sightings[,"DISENTGL"][k]          # disentangled at latest sighting
  whaletable[i,20] <- max(sightings[,"SICK"])            # ever sighted sick
  whaletable[i,21] <- sightings[,"SICK"][k]              # latest sick sighting
  whaletable[i,22] <- max(sightings[,"MEDICAL"])         # ever given medical assistance
  whaletable[i,23] <- sightings[,"MEDICAL"][k]           # latest medical sighting
  whaletable[i,24] <- max(sightings[,"DEAD"])            # ever known to be dead
  if(whaletable[i,2]==2)  whaletable[i,12] <- 0         # if a whale is male, set momhadcalf to 0
  
  if(tempid %in% vesselID){							    #if whale ID is in vessel strike data, save number of strikes
  	whaletable[i,25] <- max(vesselData[,"EventNo"][tempid==vesselID])
    if(grepl("SUPERFICIAL", vesselData[,"VesselStrikeComment"][tempid==vesselID], fixed="True")){
    	whaletable[i, 26] <- 1
    	whaletable[i,27] <- 0
    } 
    else if(grepl("SHALLOW", vesselData[,"VesselStrikeComment"][tempid==vesselID], fixed="True")){
    	whaletable[i, 26] <- 2
    	whaletable[i,27] <- 0
    } 
    else if(grepl("DEEP", vesselData[,"VesselStrikeComment"][tempid==vesselID], fixed="True")){
    	whaletable[i, 26] <- 3
    	whaletable[i, 27] <- 1
    }
    }  else  {
  	whaletable[i,25] <- 0								#record strikes as 0 if whale has not been struck
  	whaletable[i,26] <- 0								#record strike_severity as 0 if whale has not been struck
  	whaletable[i,27] <- 0								#record deepStrike as 0 if whale has not been struck
  }
  
  # we usually use 0 instead of NA so that in the regression, inapplicable whales will have 0 adjustment
  #   rather than being NA and dropped from the regression
  
  if(max(sightings[,"FENTGL"])==1)  {                    # if FENGTL, record the first year
    whaletable[i,15] <- min(sightings[sightings[,"FENTGL"]==1,"SightingYear"])
  }  else  {
    whaletable[i,15] <- NA
  }
  
  if(max(sightings[,"DISENTGL"])==1)  {                  # if DISENGTL, record the last year
    whaletable[i,18] <- max(sightings[sightings[,"DISENTGL"]==1,"SightingYear"])
  }  else  {
    whaletable[i,18] <- NA
  }
  
  if(!is.na(whaletable[i,15]) & is.na(whaletable[i,18]))  {   # assumed still entangled
    whaletable[i,19] <- finalyear - whaletable[i,15] + 1    # won't work if entangled >1 times
  }  else  {
    whaletable[i,19] <- whaletable[i,18] - whaletable[i,15] + 1
  }
  
}
colnames(whaletable) <- c("id","sex","age","firstyear","lastyear","sightings",
                          "sagEver","sagNow","feedEver","feedNow","momwcalfNow","momwcalfEver", 
                          "entglEver","entglNow","fentglyear","disentglEver","disentglNow","disentglyear",
                          "yearsentgl","sickEver","sickNow","medicalEver","medicalNow","dead","strikes", "strike_severity", "deepStrike")

write.table(whaletable, "whaletable.txt", sep="\t", row.names=F)
whaletable <- whaletable[whaletable[,"lastyear"]>=initialyear,]
allwhaletable <- whaletable                        # save a copy in case needed later

id <- whaletable[,"id"]

# delete sightings whose whale ID numbers are not in our database (e.g., whales not sighted post-1986)
whales_to_delete <- c() 
for(i in 1: numrows){
  if(!(data[i,"SightingEGNo"] %in% id)){
    whales_to_delete <- c(whales_to_delete,i)
  }
}
data <- data %>% dplyr::slice(-whales_to_delete)
numrows <- dim(data)[1]   # revised number of sightings, excluding whales not sighted post-1986



# recode sex to work better for logistic regression
female <- (whaletable[,"sex"])
female[female==2] <- -1                               # -1 = males, 1 = females
# Changed from NA to 0 so ~50 whales would not be thrown out
female[female==3] <- 0 # NA                              # unknown








####Logic for adjusting a whale's age####
# We either default to first sighting year - current year or the categorical variable  
# currentyear <- 2020
currentyear <- finalyear
temp <- whaletable
for(i in 1: nrow(whaletable)) {
  minage <- currentyear - whaletable[,"firstyear"][i]
  if (minage >= 9) {
    whaletable[,"age"][i] <- 1
  }
}
####.####



adult <- whaletable[,"age"]
adult[adult==2 | adult==3] <- -1                      # 1 = adult, -1 = calf/juvenile
adult[adult==4] <- 0   # NA                              # unknown


# do we want this?
# assign maximum (latest) age class of each whale to all sightings of that whale
tempnrows <- dim(data)[1]                                # no. of sightings currently in database
tempid <- data[,"SightingEGNo"]                          # whale ID numbers currently in database
maxageclass <- data[,"AgeClassCode"]                     # initialize
for(i in (tempnrows-1):1)  {                             # loop backwards over all sightings
  if(tempid[i]==tempid[i+1])                             # if this is the same whale
    maxageclass[i] <- maxageclass[i+1]                   #   then use the later age class
  # otherwise, this is the latest sighting of a new whale; leave its maxageclass as is
  # alternately, could just take the last age recorded, but this way is safer if later ages were misrecorded  
}

data <- cbind(data, maxageclass)

# check age class classifications
agenum <- data[,"Age"]
#agenum = factor(agenum,levels(agenum)[c(1:3,14,25,33:38,4:13, 15:24, 26:32, 39,40)])
#levels(agenum)
#table(agenum, data[,"AgeClassCode"])


momwcalfNow <- whaletable[,"momwcalfNow"]                           # 1 = (adult female) w calf
momwcalfNow[momwcalfNow==0 & female==1 & adult==1] <- -1            # -1 = adult female w/o calf
#                                                                   # 0 = (male or juvenile)

momwcalfEver <- whaletable[,"momwcalfEver"]
momwcalfEver[momwcalfEver==0 & female==1 & adult==1] <- -1

firstyear <- whaletable[,"firstyear"]
lastyear <- whaletable[,"lastyear"]
sightings <- whaletable[,"sightings"]

sagEver <- whaletable[,"sagEver"]
sagNow <- whaletable[,"sagNow"]

feedEver <- whaletable[,"feedEver"]
feedNow <- whaletable[,"feedNow"]

entglEver <- whaletable[,"entglEver"]
entglNow <- whaletable[,"entglNow"]

fentglyear <- whaletable[,"fentglyear"]

disentglEver <- whaletable[,"disentglEver"]
disentglNow <- whaletable[,"disentglNow"]

disentglyear <- whaletable[,"disentglyear"]
yearsentgl <- whaletable[,"yearsentgl"]

sick <- whaletable[,"sickEver"] | whaletable[,"sickNow"]

medical <- whaletable[,"medicalEver"] | whaletable[,"medicalNow"]

sickmed <- sick | medical

dead <- whaletable[,"dead"]
alive <- 1 - dead

sickNow <- whaletable[,"sickNow"]
sickEver <- whaletable[,"sickEver"]

strikes <- whaletable[,"strikes"]
strike_severity <- whaletable[,"strike_severity"]
table(strike_severity)
strike_severity <- factor(strike_severity, ordered=TRUE)
strikeEver <- strikes > 0
table(strike_severity)
strikeEver2 <- strikeEver + 0.0
deepStrike <- whaletable[,"deepStrike"]

# recode disentglEver to separate out never-entangled and disentangled
disentglEver[entglEver==1 & disentglEver==0] <- -1






# Try logistic regression models

# all covarites except fentglyear, disentglyear, yearsentgl, 
# and excluding feedNow, sagNOW and disentglNow, which have huge SDs
firstyear1 <- (firstyear-initialyear)
firstyear2 <- (firstyear-initialyear)^2
fit0 <- (glm(alive ~ female + entglEver + entglNow + momwcalfEver + momwcalfNow 
             + feedEver  + sagEver + disentglEver  
             + sickmed + adult + firstyear1 + firstyear2 + deepStrike, 
             family="binomial"));   summary(fit0)   
# feedNow and sagNow have huge SDs
# note alive means not confirmed dead; it doesn't mean confirmed alive

# omit sickmed, not collected for whole time period
fit1 <- (glm(alive ~ female + entglEver + entglNow + momwcalfEver + momwcalfNow 
             + feedEver  + sagEver + disentglEver  
             + adult + firstyear1 + firstyear2, 
             family="binomial"));   summary(fit1)  

# omit sagEver and feedEver since most whales have been seen doing these behaviors
fit2 <- (glm(alive ~ female + entglEver + entglNow + momwcalfEver + momwcalfNow 
             + disentglEver  
             + adult + firstyear1 + firstyear2, 
             family="binomial"));   summary(fit2)  

# Old models we have not worked on together
# fit0 <- (glm(alive ~ female + entglEver + momwcalfEver + sagNow, family="binomial"));   summary(fit0)   # AIC ...        ACTUAL(...)
# fit1 <- (glm(alive ~ female + momwcalfEver + entglEver, family="binomial"));   summary(fit1)   # AIC 194.2        ACTUAL(191.38)
# fit2 <- (glm(alive ~ adult + momwcalfNow + momwcalfEver + female + entglEver, family="binomial"));   summary(fit2)   # AIC 194.2        ACTUAL(191.38)
# fit2 <- (glm(alive ~ adult + momwcalfEver + female + entglEver + entglNow, family="binomial"));   summary(fit2)   # AIC 194.2        ACTUAL(191.38)

# Likely most robust model
fit3 <- (glm(alive ~ adult + momwcalfNow + momwcalfEver + female + entglEver + entglNow + disentglEver, family="binomial"));          summary(fit3)   # AIC 196          ACTUAL(196)

####Testing the assumptions of the model's####
# This list vector should have all the models that are being considered
modelFits <- list(fit0, fit1, fit2, fit3)
par(mfrow=c(2,2))
# This figure represents a setting of the assumptions of the model
for (i in 1:length(modelFits)) {
  #Simulate residuals in DHARMa
  res.modelFits <- simulateResiduals(modelFits[[i]], refit = F, n = 1000)
  
  #Check for uniformity
  testUniformity(res.modelFits)
  
  #Test for overdispersion
  testDispersion(res.modelFits) #this uses the dharma simulaiton output, notice the res. before models
  
  #Test for zero inflation
  testZeroInflation(res.modelFits) #this uses the dharma simulaiton output, notice the res. before models
  
  #Plot the resiuals of each variable - so we are going to plot Population
  plotResiduals(res.modelFits, asFactor = T)
  # plotResiduals(res.modelFits, asFactor = F)
  
  #None of the populations are too terrible in residuals so we could keep them
  summary(modelFits[[i]])
}
####.####



# Trying old model with all relevent variables: sex, entangle, age, w/ calf: original model below
# fit4 <- (glm(alive ~ adult + calf + female + entgl, family="binomial"));   summary(fit4)   # 174.7    ACTUAL(172.45)
fit3 <- (glm(alive ~ adult + momwcalfNow + momwcalfEver + female + entglEver + entglNow, family="binomial"));          summary(fit3)   # AIC 196          ACTUAL(196)

# Best model from the summer
fit4 <- (glm(alive ~ adult+feedEver+momwcalfNow+momwcalfEver+entglEver+entglNow, family="binomial"));   summary(fit4)   # 174.7    ACTUAL(172.45)

# Likely most robust model
# fit4 <- (glm(alive ~ adult + feedEver + momwcalfNow + momwcalfEver + entglEver + entglNow, family="binomial"));   summary(fit4)   # 174.7    ACTUAL(172.45)
# Vitor's model which is what is used in the code below 
fit3 <- (glm(alive ~ adult + momwcalfNow + momwcalfEver + female + 
               entglEver + entglNow + disentglNow, family="binomial"));          summary(fit3)   # AIC 196          ACTUAL(196)
coeffs <- coef(fit3)[2:8]
coeffs









# Process sightings database: remove associated sightings
# remove confirmed dead whales from dataset
# note: whale 1128 was already removed since all its sightings predate 1986
DEAD <- data[,"DEAD"]                         # column created by Ling for known dead whales
dead <- data[DEAD==1,"SightingEGNo"]          # list of ID numbers including duplicates
dead <- unique(dead)                          # remove duplicate ID numbers
for(i in 1:length(dead))                      # remove dead whales from sighting-level database
  data <- data[data[,"SightingEGNo"]!=dead[i],]
whaletable <- whaletable[whaletable[,"dead"]==0,]        # remove dead whales from whaletable


# remove sightings with month = 0 or day = 0  (only a few dozen of these; just delete them)
data <- data[data[,"SightingDay"]!=0,]
data <- data[data[,"SightingMonth"]!=0,]      # these should have been removed earlier (< 1986)


# read in information on potentially living whales
id <- data[,"SightingEGNo"]       # whale ID number
nwhales <- length(unique(id))        # number of sighted whales
idlist <- unique(id)                 # list of distinct whale ID numbers
numrows <- dim(data)[1]             # number of sightings in database
#   (sightings as defined by NARWC)


# convert sighting dates and times to standard format
year <- data[,"SightingYear"]                            # year of sighting
month <- data[,"SightingMonth"]                          # month of sighting
day <- data[,"SightingDay"]                              # day of sighting
time <- data[,"SightingTime"]                            # time of sighting
stdtime <- as.numeric(converttime(year,month,day,time))  # convert to standard time format
data <- cbind(data, stdtime)



# calculate waiting times in days (sighting gaps) since last sighting
# this is used only to find whales that seem to be exceeding top speed, as an error check

# first, must deal with missing times
nummissingtimes <- length(which(time==0))                # 1742 sightings have missing times
# qqnorm(time[time!=0])                 # looks reasonably normal; impute normally dist. times
imputedtimes <- rnorm(nummissingtimes, mean(time[time!=0]),sd(time[time!=0]))
timeimp <- time
timeimp[timeimp==0] <- imputedtimes


waittime <- stdtime[2:numrows] - stdtime[1:(numrows-1)]
waittime <- c(NA, waittime)           # first sighting in database has no waiting time
for(i in 2:numrows)                   # if not the same whale, then make its waittime = NA
  if (id[i]!=id[i-1])    waittime[i] <- NA
data$waittime <- round(waittime,5)   # add waiting time to sighting-level database
# check to make sure there are no negative values (confirms that database is sorted correctly)
summary(waittime)

# calculate distance traveled since last sighting (in km)
lat <- data[,"Latitude"]                                 # latitude of sighting
long <- data[,"Longitude"]                               # longitude of sighting
lat[lat==0] <- NA;    long[long==0] <- NA;        # missing data were recorded as 0; change to NA
traveldist <- rep(NA,numrows)
for(i in 2:numrows)  {
  if (id[i]==id[i-1])  {             # if this is the same whale as previous row
    lat1 <- deg2rad(lat[i-1]);    long1 <- deg2rad(long[i-1])
    lat2 <- deg2rad(lat[i]);      long2 <- deg2rad(long[i])
    traveldist[i] <- gcd.slc(long1, lat1, long2, lat2)
  }
  # if not the same whale, leave its traveldist as NA
}




data$traveldist <- round(traveldist,5)



# sightings whose traveldist/waittime combination exceeds top speed of whale
data$speed <- round(traveldist/(waittime*24),3)
maxspeed <- 9.3
toofast <- rep(0, numrows)
prevtime <- rep(0, numrows)


# merge (or delete) linked or associated sightings: within 1 month in same region code
databkp <- data                                          # save a backup copy
region <- data[,"RegionCode"]
rowstokeep <- rep(1, numrows)
waittimecutoff <- 30
traveldistcutoff <- 60
for(i in 2:numrows)  {
  if( !is.na(waittime[i]) & waittime[i]<=waittimecutoff & region[i]==region[i-1] )
    rowstokeep[i] <- 0
  
  if(!is.na(traveldist[i]) & !is.na(waittime[i]) & traveldist[i]>waittime[i]*9.3*24)
    rowstokeep[i] <- 0                   # these must be errors - exceeds top speed of whale
  toofast[i]    <- 1
  prevtime[i] <- time[i-1]
}
# this cleaning code comes after the whaletable is created, because we don't want to remove sightings
# that might contain information abotu covariates when we create whaletable
write.table(data[toofast==1,], "toofast.txt", sep="\t", row.names=F)

data <- cbind(data, prevtime)    # this could be deleted
data <- data[rowstokeep==1,]




# data <- cbind(data, maxageclass)
rm(tempid, tempnrows)                                    # clean up





# These lines are used only for selecting subsets of the data

# select by sex
sex <- data[,"GenderCode"]
# data <- data[sex=="F",]

# select by age
maxageclass <- data[,"maxageclass"]
# data <- data[maxageclass=="A",]


# clean up
rm(stdtime, waittime, traveldist, maxageclass)
# read these in fresh if they are needed later, since the number of rows may have changed





# -------------------------------- Main analysis ----------------------------- #

t.start <- Sys.time()



# update information on potentially living whales
id <- whaletable[,"id"]             # whale ID numbers
nwhales <- length(table(id))              # number of sighted whales
idlist <- rownames(table(id))                            # list of distinct whale ID numbers
stdtime <- data[,"stdtime"]                              # sighting dates

# initialize
numreps <- 100#10                                        # number of random samples to draw


# Register cores for parallel computing
cores=detectCores()
cl <- makeCluster(cores-1) # not to overload your computer
registerDoParallel(cl)

b<-nrow(whaletable)

# Parallel, create vector of P(alive for all sighted whales)
dev.off()
# Calculates probability alive using Solow method
probalive_solow <- foreach(i=1:b, .combine="c", .inorder=FALSE) %dopar% {
  tempid <- whaletable[i,"id"]                            # ID number of current whale
  tempsightings <- stdtime[id==tempid]                   # select matching sightings
  prob_i_solow <- solow(tempsightings, now=now)
  prob_i_solow
}

# Calculates probability alive using Steve's ABM method
# Unclear if this us how to accurately use this method. Should be reviewed
probalive_abm <- foreach(i=1:b, .combine="c", .inorder=FALSE) %dopar% {
  tempid <- whaletable[i,"id"]                            # ID number of current whale
  tempsightings <- stdtime[id==tempid]                   # select matching sightings
  # Changing the ext value to be either TRUE or FALSE really effects the final result.
  # F gives a distribution of estimated survival that is very similar to the Solow method
  # As of 1/9/2020 ext is set to T but I am unsure why that is the case.
  prob_i_abm <- abm38inv(tempsightings, ext=T, distance=T, now=now)
  prob_i_abm
}

# Saves a copy of the probabilities for future reference
probaliveold_solow <- probalive_solow                            
probaliveold_abm <- probalive_abm
print(sum(probaliveold_abm)) #sum of P(alive) before logistic regression


####Trying to figure out what is causing such a difference in abm and solow####
# Assumes that since abm and solow is calculated per row in a table sorted by whale 
# then the vector of probabilities from solow and abm should map onto eachother
plot(probalive_solow, probalive_abm)

# Visualize the distribution of alive in histogram
par(mfrow=c(1,2))
hist(probalive_solow, breaks = 40)
hist(probalive_abm, breaks = 40)
####.####


Sys.time() - t.start


# adult, feedEver, momwcalfNow, momwcalfEver, entglEver, entglNow
# adjust P(alive) to account for covariates
oddsalive_solow <- probalive_solow/(1-probalive_solow)
oddsalive_abm <- probalive_abm/(1-probalive_abm)


####Create tempwhaletable that associates solow and abm methods to whaletable####
tempwhaletable <- as_tibble(whaletable) %>% dplyr::mutate(solow = probalive_solow, abm = probalive_abm)
####.####



####Uses covariate information from logistic regression to adjust Solow and ABM probabilities of alive####
# Updates Solow's probability of alive given whale covariates
coeffs
for(i in 1:b)  {
  tempid <- whaletable[i,"id"]              # ID number of current whale
  currwhale <- whaletable[whaletable[,"id"]==tempid,]   # covariate information for this whale
  # Old model
  #x <- currwhale[c("age","feedEver","momwcalfNow","momwcalfEver","entglEver","entglNow")] # vector of covariates (must be in correct order)
  # Vitor's model
  x <- currwhale[c("age", "momwcalfNow","momwcalfEver", "sex", "entglEver","entglNow", "disentglNow")] # vector of covariates (must be in correct order)
  
  if(x["age"]==2 | x["age"]==3)  x["age"] <- -1     #reparametrize age
  
  if(x["age"]==4)  x["age"] <- 0
  
  if(x["sex"]==2)  x["sex"] <- -1              #reparametrize sex
  
  if(x["sex"]==3)  x["sex"] <- 0
  
  #if(x["momwcalfEver"]==0 & x["age"]==1 & x["sex"]==1)  x["momwcalfEver"] <- -1   # reparametrize calf
  if(x["momwcalfEver"]==0 & x["age"]==1)  x["momwcalfEver"] <- -1
  # other covariates are parametrized as 0-1 because they have no NA's
  adjustment <-  exp(sum(coeffs*x))
  oddsalive_solow[i] <- oddsalive_solow[i]*adjustment
}
probalive_solow <- oddsalive_solow/(1+oddsalive_solow)

# Updates ABM's probability of alive given whale covariates
coeffs
for(i in 1:b)  {
  tempid <- whaletable[i,"id"]              # ID number of current whale
  currwhale <- whaletable[whaletable[,"id"]==tempid,]   # covariate information for this whale
  # Old model
  #x <- currwhale[c("age","feedEver","momwcalfNow","momwcalfEver","entglEver","entglNow")] # vector of covariates (must be in correct order)
  # Vitor's model
  x <- currwhale[c("age", "momwcalfNow","momwcalfEver", "sex", "entglEver","entglNow", "disentglNow")] # vector of covariates (must be in correct order)
  
  if(x["age"]==2 | x["age"]==3)  x["age"] <- -1     #reparametrize age
  
  if(x["age"]==4)  x["age"] <- 0
  
  if(x["sex"]==2)  x["sex"] <- -1              #reparametrize sex
  
  if(x["sex"]==3)  x["sex"] <- 0
  
  #if(x["momwcalfEver"]==0 & x["age"]==1 & x["sex"]==1)  x["momwcalfEver"] <- -1   # reparametrize calf
  if(x["momwcalfEver"]==0 & x["age"]==1)  x["momwcalfEver"] <- -1
  # other covariates are parametrized as 0-1 because they have no NA's
  adjustment <-  exp(sum(coeffs*x))
  oddsalive_abm[i] <- oddsalive_abm[i]*adjustment
}
probalive_abm <- oddsalive_abm/(1+oddsalive_abm)
print(sum(probalive_abm)) #sum of P(alive) after logistic regression
####.####


## Histogram of old and new probability of survival
par(mfrow=c(2,2))
## We would expect a more even distribution of probability alive
## from 0-1 but until methods are fixed have different xlim scales
# hist(probaliveold_solow, breaks = 50, xlim = c(0,1))
# hist(probalive_solow, breaks = 50, xlim = c(0,1))
# hist(probaliveold_abm, breaks = 50, xlim = c(0,1))
# hist(probalive_abm, breaks = 50, xlim = c(0,1))
hist(probaliveold_solow, breaks = 50)
hist(probalive_solow, breaks = 50)
hist(probaliveold_abm, breaks = 50)
hist(probalive_abm, breaks = 50)

age <- whaletable[,"age"]
feedEver <- whaletable[,"feedEver"]
momwcalfNow <- whaletable[,"momwcalfNow"]
momwcalfEver <- whaletable[,"momwcalfEver"]
entglEver <- whaletable[,"entglEver"]
entglNow <- whaletable[,"entglNow"]
sex <- whaletable[,"sex"]


####This commented out chunk of code calculates old and new probabilities of alive####
# Old is the probability of alive just using solow or ABM
# New is adjusted probability of alive from logistic regression
# Commented out because both Solow and ABM are currently giving strange results 

# par(mfrow=c(2,3))
# 
# plot(probaliveold, probalive, type="n", xlim = c(0,1), ylim = c(0,1))
# subset1 <- whaletable
# points(probaliveold[subset1], probalive[subset1], col="red")
# abline(0,1)
# title(main="ALL WHALES")
# 
# 
# plot(probaliveold, probalive, type="n", xlim = c(0,1), ylim = c(0,1))
# subset1 <- entglEver==1 & entglNow==1
# points(probaliveold[subset1], probalive[subset1], col="red")
# subset2 <- entglEver==0
# points(probaliveold[subset2], probalive[subset2], col="blue")
# subset3 <- entglEver==1 & entglNow==0
# points(probaliveold[subset3],probalive[subset3], col = "gold")
# abline(0,1)
# legend(-.03,1.03, c("entangled now","never entangled","ever entangled"), fill=c("red","blue", "gold"), cex=.7, border=F, bty="n")
# title(main="entangled")
# 
# plot(probaliveold, probalive, type="n", xlim = c(0,1), ylim = c(0,1))
# subset1 <- sex==1
# points(probaliveold[subset1], probalive[subset1], col="red")
# subset2 <- sex==2
# points(probaliveold[subset2], probalive[subset2], col="blue")
# legend(-.03,1.03, c("female","non-female"), fill=c("red","blue"), cex=.7, border=F, bty="n")
# abline(0,1)
# title(main="sex")
# 
# plot(probaliveold, probalive, type="n", xlim = c(0,1), ylim = c(0,1))
# subset1 <- feedEver==1
# points(probaliveold[subset1], probalive[subset1], col="red")
# subset2 <- feedEver==0
# points(probaliveold[subset2], probalive[subset2], col="blue")
# abline(0,1)
# legend(-.03,1.03, c("fed","never fed"), fill=c("red","blue"), cex=.7, border=F, bty="n")
# title(main="feed")
# 
# plot(probaliveold, probalive, type="n", xlim = c(0,1), ylim = c(0,1))
# subset1 <- age==1
# points(probaliveold[subset1], probalive[subset1], col="red")
# subset2 <- age!=1
# points(probaliveold[subset2], probalive[subset2], col="blue")
# abline(0,1)
# legend(-.03,1.03, c("adult","non-adult"), fill=c("red","blue"), cex=.7, border=F, bty="n")
# title(main="age")
# 
# plot(probaliveold, probalive, type="n", xlim = c(0,1), ylim = c(0,1))
# subset1 <- momwcalfEver==1
# points(probaliveold[subset1], probalive[subset1], col="red")
# subset2 <- momwcalfEver!=1
# points(probaliveold[subset2], probalive[subset2], col="blue")
# abline(0,1)
# legend(-.03,1.03, c("momwcalfEver=1","momwcalfEver!=1"), fill=c("red","blue"), cex=.7, border=F, bty="n")
# title(main="momwcalfEver")
# 
# 
# whaletable <- cbind(whaletable, probaliveold, probalive)
# write.table(whaletable, "whaletable.txt", sep="\t", row.names=F)
# 
# 
# 
# plot(probaliveold, probalive, type="n", xlab="P(theta | X)", ylab=" P(theta | X, Y)", xlim = c(0,1), ylim = c(0,1))
# subset1 <- momwcalfEver==0 & age==1 & sex==2 & entglEver==0 & feedEver==1
# points(probaliveold[subset1], probalive[subset1], col="red", pch=16)
# subset2 <- momwcalfEver==0 & age==1 & sex==2 & entglEver==1 & feedEver==1
# points(probaliveold[subset2], probalive[subset2], col="blue", pch=16)
# subset3 <- momwcalfEver==0 & age==1 & sex==2 & entglEver==1 & feedEver==0
# points(probaliveold[subset3], probalive[subset3], col="gold", pch=16)
# subset4 <- momwcalfEver==0 & age==1 & sex==2 & entglEver==1 & feedEver==0
# points(probaliveold[subset4], probalive[subset4], col="hotpink", pch=16)
# abline(0,1, col=gray(.8))
# legend(-.03,1.03, c("M/adult/not entgl/feed","M/adult/entgl/feed","M/adult/entgl/not feed","M/adult/not entgl/not feed"), fill=c("red","blue","gold","hotpink"), cex=.7, border=F, bty="n")
# title(main="Fig. 1: Adult Males")
# 
# plot(probaliveold, probalive, type="n", xlab="P(theta | X)", ylab=" P(theta | X, Y)", xlim = c(0,1), ylim = c(0,1))
# subset1 <- momwcalfEver==0 & sex==2 & entglEver==0 & feedEver==1 & age!=1
# points(probaliveold[subset1], probalive[subset1], col="red", pch=16)
# subset2 <- momwcalfEver==0 & sex==2 & entglEver==1 & feedEver==1 & age!=1
# points(probaliveold[subset2], probalive[subset2], col="blue", pch=16)
# subset3 <- momwcalfEver==0 & sex==2 & entglEver==1 & feedEver==0 & age!=1
# points(probaliveold[subset3], probalive[subset3], col="gold", pch=16)
# subset4 <- momwcalfEver==0 & sex==2 & entglEver==1 & feedEver==0 & age!=1
# points(probaliveold[subset4], probalive[subset4], col="hotpink", pch=16)
# abline(0,1, col=gray(.8))
# legend(-.03,1.03, c("M/not entgl/feed","M/entgl/feed","M/entgl/not feed","M/not entgl/not feed"), fill=c("red","blue","gold","hotpink"), cex=.7, border=F, bty="n")
# title(main="Fig. 2: Non-adult Males")
# 
# plot(probaliveold, probalive, type="n", xlab="P(theta | X)", ylab=" P(theta | X, Y)", xlim = c(0,1), ylim = c(0,1))
# subset1 <- momwcalfEver==0 & sex==1 & entglEver==0 & feedEver==1 & age==1
# points(probaliveold[subset1], probalive[subset1], col="red", pch=16)
# subset2 <- momwcalfEver==0 & sex==1 & entglEver==1 & feedEver==1 & age==1
# points(probaliveold[subset2], probalive[subset2], col="blue", pch=16)
# subset3 <- momwcalfEver==0 & sex==1 & entglEver==1 & feedEver==0 & age==1
# points(probaliveold[subset3], probalive[subset3], col="gold", pch=16)
# subset4 <- momwcalfEver==0 & sex==1 & entglEver==1 & feedEver==0 & age==1
# points(probaliveold[subset4], probalive[subset4], col="hotpink", pch=16)
# abline(0,1, col=gray(.8))
# legend(-.03,1.03, c("F/adult/not entgl/feed","F/adult/entgl/feed","F/adult/entgl/not feed","F/adult/not entgl/not feed"), fill=c("red","blue","gold","hotpink"), cex=.7, border=F, bty="n")
# title(main="Fig. 3: Adult females (never w calf)")
# 
# 
# plot(probaliveold, probalive, type="n", xlab="P(theta | X)", ylab=" P(theta | X, Y)", xlim = c(0,1), ylim = c(0,1))
# subset1 <- momwcalfEver==1 & sex==1 & entglEver==0 & feedEver==1 & age==1
# points(probaliveold[subset1], probalive[subset1], col="red", pch=16)
# subset2 <- momwcalfEver==1 & sex==1 & entglEver==1 & feedEver==1 & age==1
# points(probaliveold[subset2], probalive[subset2], col="blue", pch=16)
# subset3 <- momwcalfEver==1 & sex==1 & entglEver==1 & feedEver==0 & age==1
# points(probaliveold[subset3], probalive[subset3], col="gold", pch=16)
# subset4 <- momwcalfEver==1 & sex==1 & entglEver==1 & feedEver==0 & age==1
# points(probaliveold[subset4], probalive[subset4], col="hotpink", pch=16)
# abline(0,1, col=gray(.8))
# legend(-.03,1.03, c("F/adult/not entgl/feed","F/adult/entgl/feed","F/adult/entgl/not feed","F/adult/not entgl/not feed"), fill=c("red","blue","gold","hotpink"), cex=.7, border=F, bty="n")
# title(main="Fig. 4: Adult females (w calf ever)")
# 
# plot(probaliveold, probalive, type="n", xlab="P(theta | X)", ylab=" P(theta | X, Y)", xlim = c(0,1), ylim = c(0,1))
# subset1 <- momwcalfEver==0 & sex==1 & entglEver==0 & feedEver==1 & age!=1
# points(probaliveold[subset1], probalive[subset1], col="red", pch=16)
# subset2 <- momwcalfEver==0 & sex==1 & entglEver==1 & feedEver==1 & age!=1
# points(probaliveold[subset2], probalive[subset2], col="blue", pch=16)
# subset3 <- momwcalfEver==0 & sex==1 & entglEver==1 & feedEver==0 & age!=1
# points(probaliveold[subset3], probalive[subset3], col="gold", pch=16)
# subset4 <- momwcalfEver==0 & sex==1 & entglEver==1 & feedEver==0 & age!=1
# points(probaliveold[subset4], probalive[subset4], col="hotpink", pch=16)
# abline(0,1, col=gray(.8))
# legend(-.03,1.03, c("F/not entgl/feed","F/entgl/feed","F/entgl/not feed","F/not entgl/not feed"), fill=c("red","blue","gold","hotpink"), cex=.7, border=F, bty="n")
# title(main="Fig. 5: Non-adult females")
# 
# 
# plot(probaliveold, probalive, type="n", xlab="P(theta | X)", ylab=" P(theta | X, Y)", xlim = c(0,1), ylim = c(0,1))
# subset1 <- entglEver==0 & feedEver==1
# points(probaliveold[subset1], probalive[subset1], col="red", pch=16)
# subset2 <- entglEver==1 & feedEver==1
# points(probaliveold[subset2], probalive[subset2], col="blue", pch=16)
# subset3 <- entglEver==1 & feedEver==0
# points(probaliveold[subset3], probalive[subset3], col="gold", pch=16)
# subset4 <- entglEver==1 & feedEver==0
# points(probaliveold[subset4], probalive[subset4], col="hotpink", pch=16)
# abline(0,1, col=gray(.8))
# legend(-.03,1.03, c("not entgl/feed","entgl/feed","entgl/not feed","not entgl/not feed"), fill=c("red","blue","gold","hotpink"), cex=.7, border=F, bty="n")
# title(main="Fig. 6: Feed and Entgl only")

####.####



cores=detectCores()
cl <- makeCluster(cores-1) # not to overload your computer
registerDoParallel(cl)

# Parallel, Random sampling using breakaway on both Solow and ABM
#    numreps, number of simulations
#   .combine = rbind, simulation outputs are stacked by row
#   .inorder = FALSE, simulations are not order dependent
results_solow = foreach(rep=1:numreps, .combine=rbind, .inorder=FALSE) %dopar% {
  
  currdata <- data                                # whales sampled dead are removed from data so we store data in temporary variable, currdata
  
  theta <- rbinom(b,1,probalive_solow)            # vector of 0 (dead) or 1 (alive) for each whale
  S <- sum(theta)                                 # number of living whales in this iteration
  
  # exclude sightings of whales for which theta = 0
  currdead <- idlist[ (1:nwhales)[theta==0]]
  for(i in 1:length(currdead))
    currdata  <- currdata[currdata[,"SightingEGNo"]!=currdead[i],]
  currid <- currdata[,"SightingEGNo"]             # update ID numbers of remaining whales
  
  # create and save 'frequency of frequencies' dataset for use by breakaway
  freqs <- table(currid)                          # number of sightings for each whale
  f <- hist(freqs, seq(min(freqs)-.5,max(freqs)+.5,1),plot=F)$counts     # frequency of frequencies
  input <- cbind(1:length(f), f)
  input <- input[input[,2]!=0, ]                                 # Removing index with 0 frequency
  input <- data.frame(index=input[,1],frequency=input[,2])    # Fixing data type for CatchAll input
  
  
  #invoke breakaway
  breakawayresults <- breakaway::breakaway(input)
  Nhat     <- breakawayresults$est
  SE       <- breakawayresults$seest
  model    <- breakawayresults$name
  tau      <- NA         # not sure what this is
  Mhat     <- Nhat - S
  phihat   <- NA
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
  
  c(tau, S, Nhat, SE, Mhat, alphahat, phihat, M, N, model)
}
####.####

cores=detectCores()
cl <- makeCluster(cores-1) # not to overload your computer
registerDoParallel(cl)
####Need to round because values too close to 0 are rounded to 0 which ruins function below####
probalive_abm <- plyr::round_any(probalive_abm, 0.1, f=ceiling)

results_abm = foreach(rep=1:numreps, .combine=rbind, .inorder=FALSE) %dopar% {
  
  currdata <- data                                # whales sampled dead are removed from data so we store data in temporary variable, currdata
  
  theta <- rbinom(b,1,probalive_abm)            # vector of 0 (dead) or 1 (alive) for each whale
  S <- sum(theta)                                 # number of living whales in this iteration
  
  # exclude sightings of whales for which theta = 0
  currdead <- idlist[ (1:nwhales)[theta==0]]
  for(i in 1:length(currdead))
    currdata  <- currdata[currdata[,"SightingEGNo"]!=currdead[i],]
  currid <- currdata[,"SightingEGNo"]             # update ID numbers of remaining whales
  
  # create and save 'frequency of frequencies' dataset for use by breakaway
  freqs <- table(currid)                          # number of sightings for each whale
  f <- hist(freqs, seq(min(freqs)-.5,max(freqs)+.5,1),plot=F)$counts     # frequency of frequencies
  input <- cbind(1:length(f), f)
  input <- input[input[,2]!=0, ]                                 # Removing index with 0 frequency
  input <- data.frame(index=input[,1],frequency=input[,2])    # Fixing data type for CatchAll input
  
  
  #invoke breakaway
  breakawayresults <- breakaway::breakaway(input)
  Nhat     <- breakawayresults$est
  SE       <- breakawayresults$seest
  model    <- breakawayresults$name
  tau      <- NA         # not sure what this is
  Mhat     <- Nhat - S
  phihat   <- NA
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
  
  c(tau, S, Nhat, SE, Mhat, alphahat, phihat, M, N, model)
}

Sys.time() - t.start
stopCluster(cl)

head(results_solow)
head(results_abm)



# Get Nhat results as numeric
Nhat_results_solow <- as.numeric(results_solow[,9])
Nhat_results_abm <- as.numeric(results_abm[,9])

# make histogram of predicted N-hat values
write.table(results_solow, "results_solow.txt", sep="\t", row.names=FALSE)
write.table(results_abm, "results_abm.txt", sep="\t", row.names=FALSE)

par(mfrow=c(1,2))
hist(Nhat_results_solow, col="darkgray", border="white", ylab="", xlab="population size",
     main="SOLOW: Estimated population size of the North Atlantic right whale", breaks = 200)#, xlim = c(250,320))

hist(Nhat_results_abm, col="darkgray", border="white", ylab="", xlab="population size",
     main="ABM: Estimated population size of the North Atlantic right whale", breaks = 200)#, xlim = c(250,320))

cat("mean of catchall predictions Nhat")
mean(na.omit(Nhat_results_solow))
mean(na.omit(Nhat_results_abm))
####.####

## Added final figure with all categories of whales
coeffs
# Variable names "age", "momwcalfNow","momwcalfEver","entglEver","entglNow", "sex"



####SOLOW: old and new prob alive####
par(mfrow=c(1,1))
plot(probaliveold_solow, probalive_solow, type="n", xlim = c(0,1), ylim = c(0,1))
subset1 <- whaletable
points(probaliveold_solow[subset1], probalive_solow[subset1], col="red", pch=16)

tmpwhaletable <- as.tibble(whaletable)

# Create subsets
subsetFemale <- tmpwhaletable$"sex"==1
subsetMale <- tmpwhaletable$"sex"==2

subsetMother <- tmpwhaletable$"sex"==1 & tmpwhaletable$"momwcalfEver"==1
subsetMotherWithCalf <- tmpwhaletable$"sex"==1 & tmpwhaletable$"momwcalfNow"==1

subsetMotherEntgl <- tmpwhaletable$"sex"==1 & tmpwhaletable$"momwcalfEver"==1 & tmpwhaletable$"entglEver"==1
subsetMotherWithCalfEntgl <- tmpwhaletable$"sex"==1 & tmpwhaletable$"momwcalfNow"==1 & tmpwhaletable$"entglEver"==1

subsetNonadult <- tmpwhaletable$"age"!=1
# subsetNotMother <- tmpwhaletable$"sex"==1 & tmpwhaletable$"momwcalfEver"==0 
# Plot points
points(probaliveold_solow[subsetFemale], probalive_solow[subsetFemale], col="blue", pch=16)
points(probaliveold_solow[subsetMale],probalive_solow[subsetMale], col = "cyan", pch=16)

points(probaliveold_solow[subsetMother], probalive_solow[subsetMother], col="darkorange", pch=16)
points(probaliveold_solow[subsetMotherWithCalf],probalive_solow[subsetMotherWithCalf], col = "gold", pch=16)

points(probaliveold_solow[subsetMotherEntgl], probalive_solow[subsetMotherEntgl], col="purple", pch=16)
points(probaliveold_solow[subsetMotherWithCalfEntgl],probalive_solow[subsetMotherWithCalfEntgl], col = "palevioletred1", pch=16)

# points(probaliveold_solow[subsetNonadult], probalive_solow[subsetNonadult], col="black", pch=16)

legend(-.03,1.00, c("Female", "Male", "Mother entangled", "Mother w/ calf", "Mother", "Mother w/ calf entangled", "Other whales"), fill=c("blue","cyan", "darkorange", "gold", "purple", "palevioletred1", "Red"), cex=.7, border=F, bty="n")

abline(0,1)
title(main="Figure 1: Effect of adjusting for covariates")
####.####



####ABM: old and new prob alive####
plot(probaliveold_abm, probalive_abm, type="n", xlim = c(0,1), ylim = c(0,1))
subset1 <- whaletable
#points(probaliveold_abm[subset1], probalive_abm[subset1], col="red", pch=1)

tmpwhaletable <- as.tibble(whaletable)

# Create subsets
subsetFemale <- tmpwhaletable$"sex"==1 & tmpwhaletable$"momwcalfEver"==0 & tmpwhaletable$"momwcalfNow"==0 & tmpwhaletable$"entglEver"==0
subsetMale <- tmpwhaletable$"sex"==2

subsetMother <- tmpwhaletable$"sex"==1 & tmpwhaletable$"momwcalfEver"==1 & tmpwhaletable$"entglEver"==0
subsetMotherWithCalf <- tmpwhaletable$"sex"==1 & tmpwhaletable$"momwcalfNow"==1 & tmpwhaletable$"entglEver"==0

subsetMotherEntgl <- tmpwhaletable$"sex"==1 & tmpwhaletable$"momwcalfEver"==1 & tmpwhaletable$"entglEver"==1 
subsetMotherWithCalfEntgl <- tmpwhaletable$"sex"==1 & tmpwhaletable$"momwcalfNow"==1 & tmpwhaletable$"entglEver"==1

subsetNonadult <- tmpwhaletable$"age"!=1
# subsetNotMother <- tmpwhaletable$"sex"==1 & tmpwhaletable$"momwcalfEver"==0 
# Plot points
points(probaliveold_abm[subsetFemale], probalive_abm[subsetFemale], col="blue", pch=2, cex=.5)
points(probaliveold_abm[subsetMale],probalive_abm[subsetMale], col = "cyan", pch=1, cex=.5)

points(probaliveold_abm[subsetMother], probalive_abm[subsetMother], col="darkorange", pch=2, cex=.5)
points(probaliveold_abm[subsetMotherWithCalf],probalive_abm[subsetMotherWithCalf], col = "gold", pch=2, cex=.5)

points(probaliveold_abm[subsetMotherEntgl], probalive_abm[subsetMotherEntgl], col="purple", pch=2, cex=.5)
points(probaliveold_abm[subsetMotherWithCalfEntgl],probalive_abm[subsetMotherWithCalfEntgl], col = "palevioletred1", pch=2, cex=.5)

# points(probaliveold_abm[subsetNonadult], probalive_abm[subsetNonadult], col="black", pch=16)

legend(-.03,1.00, c("Female", "Male", "Mother entangled", "Mother w/ calf", "Mother", "Mother w/ calf entangled", "Other whales"), fill=c("blue","cyan", "darkorange", "gold", "purple", "palevioletred1", "Red"), cex=.7, border=F, bty="n")

abline(0,1)
title(main="Figure 1: Effect of adjusting for covariates")
####.####