####What code does:####
# Code is designed to help facilitate determining which of 216 behaviores count 
# as one of the 7 categories of interest for the regression. For example, feeding 
# has more than one behavior coded that should count as feeding.
#
# Objects and what they are:
# behaviors.counts.wide -->
# behaviores.key --> read in excel of behavior key with their shorthand name, description, and comments
# behaviorTallyTable --> tally of how many times a specific coded behavior appears in the large dataset
# database --> read in observation level version of the excel dataset provided by right whale consortium
####.####

# Set your working directory with the copies of the dataset --> setwd()
#Clears out your memory (will erase EVERYTHING)
rm(list=ls())

#setwd("~/Documents/Research/Right whales/Database 2020")
setwd("~/Desktop/2020 Summer Research")

####Load packages####
library(tidyverse)
library(readxl)
library(stringr)

####Load in the whale datasets and initialize the vectors that will be used to clean dataset####
####Switch from new to old DB####
# Loaded in the updated whale dataset
behaviors.counts.wide <- read_excel("Right Whale Catalog- Feb 3 2020- for Wang and Roberts.xlsx")
# Loaded in old whale dataset
#behaviors.counts.wide <- read_excel("Right Whale Catalog 6_14_2011 (raw).xls")
####.####
behaviors.counts.wide <-behaviors.counts.wide[order(behaviors.counts.wide[,"SightingEGNo"],behaviors.counts.wide[,"SightingYear"],behaviors.counts.wide[,"SightingMonth"],behaviors.counts.wide[,"SightingDay"], behaviors.counts.wide[,"SightingTime"]),]

####Table that countes the number of times a behavior is recorded####
behaviorTallyTable <- unlist(strsplit(behaviors.counts.wide$Behaviors, ",")) 
behaviorTallyTable <- as_tibble(table(behaviorTallyTable))
behaviorTallyTable <- droplevels(behaviorTallyTable)

# Old dataset seemed to set each behavior for each observation and not keep it constant 
# Creating placeholder variables for each observation
behaviors.counts.wide$SAG          <- 0
behaviors.counts.wide$FEED         <- 0
behaviors.counts.wide$MOMWCALF     <- 0
behaviors.counts.wide$MOMHADCALF   <- 0
behaviors.counts.wide$FENTGL       <- 0
behaviors.counts.wide$NOWENTGL     <- 0
behaviors.counts.wide$EVERENTGL    <- 0
behaviors.counts.wide$DISENTGL     <- 0
behaviors.counts.wide$DEAD         <- 0
behaviors.counts.wide$SICK         <- 0
behaviors.counts.wide$MEDICAL      <- 0

# Loaded in the dataset that defines behaviores which we will use to specify which key words 
# mark a behavior as true
behaviors.key <- read_excel("Behaviors 2020-02-12.xlsx")

####Initalizing vectors that will contain behavior codes####
vectorFEED            <- as.character() 
vectorSAG             <- as.character()
vectorMOMWCALF        <- as.character()   #mom with calf
vectorMOMHADCALF      <- as.character()   #mom had calf
vectorFENTGL          <- as.character() 
vectorNOWENTGL        <- as.character()   #currently entangled
vectorEVERENTGL       <- as.character()   #ever entagled
vectorDISENTGL        <- as.character()
vectorDEAD            <- as.character()
vectorSICK            <- as.character()
vectorMEDICAL         <- as.character()
####.####

####Parse through definition and save them as vectors####
# This found all the occurances of feed in the keys behavior description and appended it to a vector
# that we use to 

####Comment: Parse through####
## Obtaining feed behavior vector
vectorTemp <- behaviors.key$Comment[grep("feed", behaviors.key$Comment, ignore.case=T)]
for (row in 1:nrow(behaviors.key)) {
    if(behaviors.key$Comment[row] %in% vectorTemp) {
        vectorFEED[length(vectorFEED)+1] <- as.character(behaviors.key$ShortName[row])
    }
}
vectorTemp <- NA

## Obtaining sag behavior vector
vectorTemp <- behaviors.key$Comment[grep("sag", behaviors.key$Comment, ignore.case=T)]
for (row in 1:nrow(behaviors.key)) {
    if(behaviors.key$Comment[row] %in% vectorTemp) {
        vectorSAG[length(vectorSAG)+1] <- as.character(behaviors.key$ShortName[row])
    }
}
vectorTemp <- NA

## Obtaining calf behavior vector for mom with calf
vectorTemp <- behaviors.key$Comment[grep("calf", behaviors.key$Comment, ignore.case=T)]
for (row in 1:nrow(behaviors.key)) {
    if(behaviors.key$Comment[row] %in% vectorTemp) {
        vectorMOMWCALF[length(vectorMOMWCALF)+1] <- as.character(behaviors.key$ShortName[row])
    }
}
vectorTemp <- NA

## Obtaining calf behavior vector for mom had calf
vectorTemp <- behaviors.key$Comment[grep("calf", behaviors.key$Comment, ignore.case=T)]
for (row in 1:nrow(behaviors.key)) {
    if(behaviors.key$Comment[row] %in% vectorTemp) {
        vectorMOMHADCALF[length(vectorMOMHADCALF)+1] <- as.character(behaviors.key$ShortName[row])
    }
}
vectorTemp <- NA

## Obtaining fentgl behavior vector
vectorTemp <- behaviors.key$Comment[grep("First Entangled", behaviors.key$Comment, ignore.case=T)]
for (row in 1:nrow(behaviors.key)) {
    if(behaviors.key$Comment[row] %in% vectorTemp) {
        vectorFENTGL[length(vectorFENTGL)+1] <- as.character(behaviors.key$ShortName[row])
    }
}
vectorTemp <- NA

## Obtaining entgl behavior vector for now entangled
vectorTemp <- behaviors.key$Comment[grep("Entangled", behaviors.key$Comment, ignore.case=T)]
for (row in 1:nrow(behaviors.key)) {
    if(behaviors.key$Comment[row] %in% vectorTemp) {
        vectorNOWENTGL[length(vectorNOWENTGL)+1] <- as.character(behaviors.key$ShortName[row])
    }
}
vectorTemp <- NA

## Obtaining entgl behavior vector for ever entangled
vectorTemp <- behaviors.key$Comment[grep("Entangled", behaviors.key$Comment, ignore.case=T)]
for (row in 1:nrow(behaviors.key)) {
    if(behaviors.key$Comment[row] %in% vectorTemp) {
        vectorEVERENTGL[length(vectorEVERENTGL)+1] <- as.character(behaviors.key$ShortName[row])
    }
}
vectorTemp <- NA

## Obtaining disentgl behavior vector for disentangled
vectorTemp <- behaviors.key$Comment[grep("disentangle", behaviors.key$Comment, ignore.case=T)]
for (row in 1:nrow(behaviors.key)) {
    if(behaviors.key$Comment[row] %in% vectorTemp) {
        vectorDISENTGL[length(vectorDISENTGL)+1] <- as.character(behaviors.key$ShortName[row])
    }
}
vectorTemp <- NA

## Obtaining dead behavior vector
vectorTemp <- behaviors.key$Comment[grep("dead", behaviors.key$Comment, ignore.case=T)]
for (row in 1:nrow(behaviors.key)) {
    if(behaviors.key$Comment[row] %in% vectorTemp) {
       vectorDEAD[length(vectorDEAD)+1] <- as.character(behaviors.key$ShortName[row])
    }
}
vectorTemp <- NA

## Obtaining sick behavior vector for sick
vectorTemp <- behaviors.key$Comment[grep("sick", behaviors.key$Comment, ignore.case=T)]
for (row in 1:nrow(behaviors.key)) {
  if(behaviors.key$Comment[row] %in% vectorTemp) {
    vectorSICK[length(vectorSICK)+1] <- as.character(behaviors.key$ShortName[row])
  }
}
vectorTemp <- NA

## Obtaining medical behavior vector
vectorTemp <- behaviors.key$Comment[grep("medical", behaviors.key$Comment, ignore.case=T)]
for (row in 1:nrow(behaviors.key)) {
  if(behaviors.key$Comment[row] %in% vectorTemp) {
    vectorMEDICAL[length(vectorMEDICAL)+1] <- as.character(behaviors.key$ShortName[row])
  }
}
vectorTemp <- NA


####Description: Parse through####
## Obtaining feed behavior vector
vectorTemp <- behaviors.key$Description[grep("feed", behaviors.key$Description, ignore.case=T)]
for (row in 1:nrow(behaviors.key)) {
  if(behaviors.key$Description[row] %in% vectorTemp) {
    vectorFEED[length(vectorFEED)+1] <- as.character(behaviors.key$ShortName[row])
  }
}
vectorTemp <- NA

## Obtaining sag behavior vector
vectorTemp <- behaviors.key$Description[grep("sag", behaviors.key$Description, ignore.case=T)]
for (row in 1:nrow(behaviors.key)) {
  if(behaviors.key$Description[row] %in% vectorTemp) {
    vectorSAG[length(vectorSAG)+1] <- as.character(behaviors.key$ShortName[row])
  }
}
vectorTemp <- NA

## Obtaining calf behavior vector for mom with calf
vectorTemp <- behaviors.key$Description[grep("calf", behaviors.key$Description, ignore.case=T)]
for (row in 1:nrow(behaviors.key)) {
  if(behaviors.key$Description[row] %in% vectorTemp) {
    vectorMOMWCALF[length(vectorMOMWCALF)+1] <- as.character(behaviors.key$ShortName[row])
  }
}
vectorTemp <- NA

## Obtaining calf behavior vector for mom had calf
vectorTemp <- behaviors.key$Description[grep("calf", behaviors.key$Description, ignore.case=T)]
for (row in 1:nrow(behaviors.key)) {
  if(behaviors.key$Description[row] %in% vectorTemp) {
    vectorMOMHADCALF[length(vectorMOMHADCALF)+1] <- as.character(behaviors.key$ShortName[row])
  }
}
vectorTemp <- NA

## Obtaining fentgl behavior vector
vectorTemp <- behaviors.key$Description[grep("First Entangled", behaviors.key$Description, ignore.case=T)]
for (row in 1:nrow(behaviors.key)) {
  if(behaviors.key$Description[row] %in% vectorTemp) {
    vectorFENTGL[length(vectorFENTGL)+1] <- as.character(behaviors.key$ShortName[row])
  }
}
vectorTemp <- NA

## Obtaining entgl behavior vector for now entangled
vectorTemp <- behaviors.key$Description[grep("Entangled", behaviors.key$Description, ignore.case=T)]
for (row in 1:nrow(behaviors.key)) {
  if(behaviors.key$Description[row] %in% vectorTemp) {
    vectorNOWENTGL[length(vectorNOWENTGL)+1] <- as.character(behaviors.key$ShortName[row])
  }
}
vectorTemp <- NA

## Obtaining entgl behavior vector for ever entangled
vectorTemp <- behaviors.key$Description[grep("Entangled", behaviors.key$Description, ignore.case=T)]
for (row in 1:nrow(behaviors.key)) {
  if(behaviors.key$Description[row] %in% vectorTemp) {
    vectorEVERENTGL[length(vectorEVERENTGL)+1] <- as.character(behaviors.key$ShortName[row])
  }
}
vectorTemp <- NA

## Obtaining disentgl behavior vector for disentangled
vectorTemp <- behaviors.key$Description[grep("disentangle", behaviors.key$Description, ignore.case=T)]
for (row in 1:nrow(behaviors.key)) {
  if(behaviors.key$Description[row] %in% vectorTemp) {
    vectorDISENTGL[length(vectorDISENTGL)+1] <- as.character(behaviors.key$ShortName[row])
  }
}
vectorTemp <- NA

## Obtaining dead behavior vector
vectorTemp <- behaviors.key$Description[grep("dead", behaviors.key$Description, ignore.case=T)]
for (row in 1:nrow(behaviors.key)) {
  if(behaviors.key$Description[row] %in% vectorTemp) {
    vectorDEAD[length(vectorDEAD)+1] <- as.character(behaviors.key$ShortName[row])
  }
}
vectorTemp <- NA

## Obtaining sick behavior vector for sick
vectorTemp <- behaviors.key$Description[grep("sick", behaviors.key$Description, ignore.case=T)]
for (row in 1:nrow(behaviors.key)) {
  if(behaviors.key$Description[row] %in% vectorTemp) {
    vectorSICK[length(vectorSICK)+1] <- as.character(behaviors.key$ShortName[row])
  }
}
vectorTemp <- NA

## Obtaining medical behavior vector
vectorTemp <- behaviors.key$Description[grep("medical", behaviors.key$Description, ignore.case=T)]
for (row in 1:nrow(behaviors.key)) {
  if(behaviors.key$Description[row] %in% vectorTemp) {
    vectorMEDICAL[length(vectorMEDICAL)+1] <- as.character(behaviors.key$ShortName[row])
  }
}
vectorTemp <- NA

####Removes duplicate values####
vectorFEED           <- unique(vectorFEED)
vectorSAG            <- unique(vectorSAG)
vectorMOMWCALF       <- unique(vectorMOMWCALF)
vectorMOMHADCALF     <- unique(vectorMOMHADCALF)
vectorFENTGL         <- unique(vectorFENTGL)
vectorNOWENTGL       <- unique(vectorNOWENTGL)
vectorEVERENTGL      <- unique(vectorEVERENTGL)
vectorDISENTGL       <- unique(vectorDISENTGL)
vectorDEAD           <- unique(vectorDEAD)
vectorSICK           <- unique(vectorSICK)
vectorMEDICAL        <- unique(vectorMEDICAL)
####.####

####Define the vectors with values that are relavant by manually looking through the automated results####
####Determine which values are relavant when looking at the descriptions of the behavior####
vectorMOMWCALF      <- c("CRDLE", "LOST CALF", "W/CALF UNPH", "W/CALF")  #add in w/calf manually
vectorMOMHADCALF    <- c("CRDLE", "LOST CALF", "W/CALF UNPH", "W/CALF")  #add in w/calf manually
vectorSAG           <- c("SAG", "APPR", "FCL")
vectorFEED          <- c("NOD", "SKM FD", "ECH", "LEAD", "SIDE FD", "FEED", "CO FD", "W/BSK SHRK", "SUB FD")
vectorFENTGL        <- c("FENTGL")
vectorNOWENTGL      <- c("PRT DSENTGL", "DSENTGL ATT", "ENTGL", "FRST ENTGL", "NOT FL", "FL")
####A whale that has ever been disentangled must have been entangled before
vectorEVERENTGL     <- c("DSENTGL", "PRT DSENTGL", "DSENTGL ATT", "ENTGL", "FRST ENTGL", "NOT FL", "FL", "LN GONE")
vectorDISENTGL      <- c("DSENTGL")
vectorDEAD          <- c("DEAD ON BEACH", "FLTG DEAD", "FRST DEAD", "RETRVD", "MORT DATA")
vectorSICK          <- c("SICK")
vectorMEDICAL       <- c("MEDICAL")
####.####

# Make all NA behavior values into an empty string
behaviors.counts.wide[is.na(behaviors.counts.wide)] <- "."

wcalfchecker <- 0
whaleID <- behaviors.counts.wide[1,"SightingEGNo"]
####Determine which behaviors an individual has exhibited based on Behavior column from raw data#####
for (row in 1:nrow(behaviors.counts.wide)) {
  if (whaleID != behaviors.counts.wide[row,"SightingEGNo"])
  {
    wcalfchecker <- 0
    whaleID <- behaviors.counts.wide[row,"SightingEGNo"]
  }
  if(wcalfchecker == 1){
    behaviors.counts.wide[row,"MOMHADCALF"] <- 1
  }
  
  # Creates a vector of an individuals behaviors split up by a comma
  thisRowBehavior <- unlist(strsplit(as.character(behaviors.counts.wide[row,"Behaviors"]), ", ", fixed = TRUE))
  if (!(length(thisRowBehavior)==1 && thisRowBehavior[1]==".")) {      
      for (i in 1:length(thisRowBehavior)) {
          # Uses vectorFEED and the Behaviores column from the uncleaned dataset (behaviors.counts.wide) to
          # determine if individuals have fed and does the same for all other behaviors
          if (thisRowBehavior[i] %in% vectorFEED) {
              behaviors.counts.wide[row,"FEED"] <- 1
          }
            
          if (thisRowBehavior[i] %in% vectorMOMWCALF) {
              behaviors.counts.wide[row,"MOMWCALF"] <- 1
              wcalfchecker <- 1
          }
          
          if (thisRowBehavior[i] %in% vectorMOMHADCALF) {
              behaviors.counts.wide[row,"MOMHADCALF"] <- 1
          }
      
          if (thisRowBehavior[i] %in% vectorDEAD) {
              behaviors.counts.wide[row,"DEAD"] <- 1
          }
        
          if (thisRowBehavior[i] %in% vectorDISENTGL) {
              behaviors.counts.wide[row,"DISENTGL"] <- 1
          }
          
          if (thisRowBehavior[i] %in% vectorNOWENTGL) {
              behaviors.counts.wide[row,"NOWENTGL"] <- 1
          }
          
          if (thisRowBehavior[i] %in% vectorEVERENTGL) {
              behaviors.counts.wide[row,"EVERENTGL"] <- 1
          }
          
          if (thisRowBehavior[i] %in% vectorFENTGL) {
              behaviors.counts.wide[row,"FENTGL"] <- 1
          }
          
          if (thisRowBehavior[i] %in% vectorSAG) {
              behaviors.counts.wide[row,"SAG"] <- 1
          }

          if (thisRowBehavior[i] %in% vectorSICK) {
              behaviors.counts.wide[row,"SICK"] <- 1
          }
          
          if (thisRowBehavior[i] %in% vectorMEDICAL) {
              behaviors.counts.wide[row,"MEDICAL"] <- 1
          }  
       }
    }
}
####.####

####Print out key for behaviors from vectors####
vectorFEED     
vectorSAG     
vectorMOMWCALF
vectorMOMHADCALF     
vectorFENTGL   
vectorNOWENTGL    
vectorEVERENTGL
vectorDISENTGL 
vectorDEAD     
vectorSICK
vectorMEDICAL
####.####

####Removing 'Behaviors' column and saving the reformatted dataset as an object file####
behaviors.counts.wide$Behaviors <- NULL
behaviors.counts.wide$SightingLetter <- NULL
behaviors.counts.wide$SightingId <- NULL

####Sort database by SightingEGNo, year, month, and day####
database <- behaviors.counts.wide %>% arrange(SightingEGNo, SightingYear, SightingMonth, SightingDay, SightingTime)
#database <- behaviors.counts.wide[order(behaviors.counts.wide[,"SightingEGNo"],behaviors.counts.wide[,"SightingYear"],behaviors.counts.wide[,"SightingMonth"],behaviors.counts.wide[,"SightingDay"], behaviors.counts.wide[,"SightingTime"]),]

write.table(database, file = "databasetemp.txt", append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)

save(database, file = "databasetemp.RData")
####Saving a sighting level temporary database with the following distinctions####
####Whale ever having a calf vs currently having a calf####
####Whale ever being entangled vs currently being entangled####
####.####