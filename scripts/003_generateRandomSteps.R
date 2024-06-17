startTime = Sys.time()

## Generate random steps from Sierra bighorn movement data
## for use in step selection functions
## Takes about 8 minutes

#### Load packages ----
library(lubridate)
library(ggplot2)
library(RStoolbox)
library(tidyverse)
library(sf)
library(rgdal)
library(terra)
library(amt)
library(purrr)
source("scripts/000_samplingParameters.R") # Seed is set inside samplingParameters.R


#### Import data ----
## Import Herd Unit shapefiles
HU_raw = readOGR("../SNBS_HerdUnits/Cathedral/",
                 "snbs_herdunits_Cathedral_BlkDvd",
                 verbose = FALSE)

# Bring in sheepwise migration date windows
dates_df = read_csv("data/mig_PatternsCompiled.csv")

## Movement data
# Import movement data
sheep_imported = read_csv("../data_ts/SNBS_Locs/SNBS_Herd_Sex_20231115_KAClean.csv") %>% 
  ## Remove the faulty Maryland data point
  ## Fix missing herd units and sex for S182, S508, and S213
  filter(Longitude < -80) %>%
  mutate(Herd_Rels = case_when(Animal_ID == "S182" ~ "Wh",
                               Animal_ID == "S508" ~ "Wh",
                               Animal_ID == "S213" ~ "Bx",
                               TRUE ~ Herd_Rels)) %>%
  mutate(Sheep_Sex = case_when(Animal_ID == "S182" ~ "F",
                               Animal_ID == "S508" ~ "M",
                               Animal_ID == "S213" ~ "F",
                               TRUE ~ Sheep_Sex))
## Convert to simple features
allsheep = st_as_sf(sheep_imported,
                    coords = c("Longitude", "Latitude"),
                    crs = 4326) %>%
  st_transform(32611) %>%
  mutate(UTME = st_coordinates(.)[,1],
         UTMN = st_coordinates(.)[,2]) %>%
  select(-c(UTM_E,UTM_N))


#### Define extraction function ----
#### This function takes a year ID and herd unit ID
#### selects individuals from the migration phenology dataset
#### resamples to the appropriate temporal frequency
#### generates a movement track
#### and returns a tibble with a herd unit - year combination's 
#### set of Individual ID's, sex, movement class, tracks and 
#### steplength + turning angle data.
#### Note that the herd unit/year breakdown is required because 
#### phenology data vary by year
## sheep_data is gps relocation data (sf object)
## x is a dash-separated character vector composed of myHU and myYear.
## myYear is a numeric of year
## myHU is a character indicator of herd unit, length 2 and capitalization Xy
## windowStyle should be "Fixed" to indicate a desired fixed width window of width windowWidth
## otherwise it will take the start and endpoints of the migratory period, which 
## can be as brief as 1 day in direct migrant cases.
randomStepsPrep = function(x = "Wh-2013",
                           sheep_data = allsheep,
                           minInt = minHr,
                           windowStyle = windStyl,
                           windowWidth = windWide){
  ## Getting started. Split up herd unit/year information
  huYr = x
  myHU_Year_split = str_split(huYr,"-")[[1]]
  myHU = myHU_Year_split[1]
  myYear = myHU_Year_split[2]
  message(paste0("Year: ", myYear))
  message(paste0("Herd Unit: ", myHU))
  
  # Reference date for snow phenology 
  # (this is the true DOY when "DOY" = 1 in snowDescriber)
  startDate = as.Date(paste0((as.numeric(myYear)-1),"-08-15"))
  
  ## Select migration phenology data
  ## If no migrants were detected in a given herd unit and year, 
  ## residents are excluded because a mean herd-unit/year migration
  ## date cannot be calculated.
  thisYrMigrantDates = dates_df %>%
    filter(springYr == myYear,
           HerdUnit == myHU) %>%
    filter(!is.na(uphillMvmt)) %>%
    dplyr::rename("Year" = springYr,
                  "midDate" = uphillMvmt,
                  "midDOY" = upDOY) %>%
    mutate(startDOY = midDOY-(0.5*uphillRate),
           stopDOY = midDOY+(0.5*uphillRate),
           startDate = as.Date(paste0(Year,"-",startDOY),format="%Y-%j"),
           stopDate = as.Date(paste0(Year,"-",stopDOY),format="%Y-%j"),
           fixedWindowStartDate = midDate - days(windowWidth/2),
           fixedWindowStopDate = midDate + days(windowWidth/2)) %>% 
    dplyr::select(Year, HerdUnit, ID, mvmtClass,
                  midDate, midDOY,
                  startDOY, stopDOY,
                  startDate, stopDate,
                  fixedWindowStartDate,
                  fixedWindowStopDate)
  
  ## If there are migrants in this herd unit and year, 
  # find appropriate dates for residents as well
  if(nrow(thisYrMigrantDates) > 0){
    # Find midpoint of herdunit migration for the year
    midMigDate = mean(thisYrMigrantDates$midDOY)
    # Apply to residents and non-upward-moving-individuals as well
    thisYrResidentDates = dates_df %>%
      filter(is.na(uphillMvmt)) %>%
      filter(springYr == myYear,
             HerdUnit == myHU) %>%
      dplyr::rename("Year" = springYr) %>%
      mutate(midDate = as.Date(paste0(myYear,"-",round(midMigDate)),
                               format = "%Y-%j"),
             midDOY = as.numeric(strftime(midDate, format = "%j"))) %>%
      mutate(startDOY = NA,
             stopDOY = NA,
             startDate = NA,
             stopDate = NA,
             fixedWindowStartDate = midDate - days(windowWidth/2),
             fixedWindowStopDate = midDate + days(windowWidth/2)) %>% 
      dplyr::select(Year, HerdUnit, ID, mvmtClass,
                    midDate, midDOY,
                    startDOY, stopDOY,
                    startDate, stopDate,
                    fixedWindowStartDate,
                    fixedWindowStopDate)  
  }else{
    ## If there are not migrants in this herd unit and year, do not try and 
    ## estimate a "migration season"
    thisYrResidentDates = NULL  
  }
  
  ## Combine the migrants and nonmigrants into one table
  thisYrDates = rbind(thisYrMigrantDates,
                      thisYrResidentDates)
  
  ## Only run the function if there are selected individuals in myHU in myYear.
  ## Technically there should be no combinations of HU-Year with no individuals
  ## because we're generating the HU-Year combinations from known animal data
  ## but it's worth checking anyway.
  if(nrow(thisYrDates) == 0){
    message("No individuals were retained in this herd unit/year.")
    message("nrow(thisYrDates == 0)")
    return(NULL)
  }
  
  ## Identify the minimum and maximum dates to be sampled for this HU/Year
  if(windowStyle == "Fixed"){
    thisYrMinDate = min(c(
      thisYrDates$fixedWindowStartDate,thisYrDates$fixedWindowStopDate),
      na.rm=T)
    thisYrMaxDate = max(c(
      thisYrDates$fixedWindowStartDate,thisYrDates$fixedWindowStopDate),
      na.rm=T)
  }else{
    thisYrMinDate = min(c(
      thisYrDates$startDate,thisYrDates$stopDate),
      na.rm=T)
    thisYrMaxDate = max(c(
      thisYrDates$startDate,thisYrDates$stopDate),
      na.rm=T)
  }
  
  ## Identify individuals in herd unit with timed uphill movement behaviors
  mySheep = thisYrDates %>%
    pull(ID)
  ## Filter sheep_data to Herd Unit
  sheep_sf = sheep_data %>%
    filter(Herd_Rels == myHU)
  ## A function to select just the dates of interest for HU-Year combinations
  gpsRetainer = function(myInd){
    goodDatesDF = thisYrDates %>%
      filter(ID == myInd,
             Year == myYear)
    if(nrow(goodDatesDF) == 0){
      return(NULL)
    }else{
      if(windowStyle == "Fixed"){
        goodDates = unlist(lapply(1:nrow(goodDatesDF),
                                  FUN = function(i){
                                    seq.Date(
                                      from=goodDatesDF[i,]$fixedWindowStartDate,
                                      to=goodDatesDF[i,]$fixedWindowStopDate,
                                      by = "1 day") %>%
                                      as.character()
                                  }))
      }else{
        goodDates = unlist(lapply(1:nrow(goodDatesDF),
                                  FUN = function(i){
                                    seq.Date(
                                      from=goodDatesDF[i,]$startDate,
                                      to=goodDatesDF[i,]$stopDate,
                                      by = "1 day") %>%
                                      as.character()
                                  }))
      }
      
      # Identify all dates of interest for this individual
      thisSheep = sheep_sf %>%
        dplyr::mutate(Date = as.character(as.Date(lubridate::mdy_hms(DateTimePST)))) %>%
        filter(Animal_ID == myInd)
      
      outDat = thisSheep %>% 
        filter(Date %in% goodDates)
      return(outDat)
    }
  }
  ## Subset all relevant individual sheep movement data from the mySheep dataset
  sheepRetainedList = lapply(mySheep,
                             gpsRetainer)
  ## Recombine all selected movement data into one data.frame()
  sheepRetained = do.call("rbind",sheepRetainedList)
  
  ## If we didn't retain any individuals in this HU-Year, return nothing
  ## (Again, shouldn't happen, but just to be sure)
  if(nrow(sheepRetained) == 0){
    message("No individuals were retained in this herd unit/year.")
    message("nrow(sheepRetained == 0) after selecting season-specific movement data.")
    return(NULL)
  }
  
  ## > Prep sheep data for amt workflow ----
  # Create a new data frame using the final selected bighorn dataset
  sheep_raw = sheepRetained %>%
    st_drop_geometry()
  
  ## Wrangle time
  # Convert to lubridate-friendly
  sheep_raw$ymdhms = mdy_hms(sheep_raw$DateTimePST)
  firstLoc = min(as.Date(sheep_raw$ymdhms))
  lastLoc = max(as.Date(sheep_raw$ymdhms))
  ## Wrangle names
  sheep_raw$x = sheep_raw$UTME
  sheep_raw$y = sheep_raw$UTMN
  sheep_raw$t = sheep_raw$ymdhms
  
  ## Format for multiple individuals (list-style)
  # Generate list of data for multiple individual analysis
  dat_all = sheep_raw %>% 
    mutate(DOY = format(ymdhms, "%j"),
           yrMin1 = as.numeric(Year) + 1,
           springYr = case_when(DOY > 304 ~ yrMin1,
                                TRUE ~ as.numeric(Year))) %>%
    mutate(IDyr = paste0(Animal_ID,"_",springYr)) %>%
    nest(-IDyr) %>%
    separate(IDyr,into=c("ID","Year"),remove = FALSE)
  
  ## Add individual-level parameters to list
  # Extract sex from raw data and assign to list
  Sex = sheep_raw %>%
    dplyr::select(Animal_ID,Sheep_Sex) %>%
    dplyr::group_by(Animal_ID,Sheep_Sex) %>%
    dplyr::summarize(IDsex = unique(Sheep_Sex))
  dat_all$Sex = Sex$IDsex[match(dat_all$ID, Sex$Animal_ID)]
  
  # Extract herd unit from raw data and assign to list
  HU = sheep_raw %>%
    group_by(Animal_ID) %>%
    dplyr::summarize(IDHU = unique(Herd_Rels))
  dat_all$HU = HU$IDHU[match(dat_all$ID, HU$Animal_ID)]
  # Extract migration class and assign to list
  mvmtClass = thisYrDates  %>%
    mutate(IDyr = paste0(ID,"_",Year))
  dat_all$MvmtClass = mvmtClass$mvmtClass[match(dat_all$IDyr, mvmtClass$IDyr)]
  
  ## > Create track ----
  # Convert to a track object
  # Retain only cases where there are data spanning the preceding
  # and succeeding 1/2 window width
  tr1 = dat_all %>%
    mutate(trk = lapply(data,
                        function(d){
                          d$DOY = as.numeric(d$DOY)
                          indID = d$Animal_ID[1]
                          minWindowDate = thisYrDates %>%
                            filter(ID == indID) %>%
                            pull(midDOY) - (windowWidth/2)
                          maxWindowDate = thisYrDates %>%
                            filter(ID == indID) %>%
                            pull(midDOY) + (windowWidth/2)
                          if(min(d$DOY) <= minWindowDate &
                             max(d$DOY) >= maxWindowDate){
                            amt::make_track(d,x,y,t,
                                            crs = st_crs(sheep_data))
                          }else(return(NULL))
                        }))
  ## Remove null tracks
  tr1 = tr1[(!unlist(lapply(tr1$trk,
                            is.null))),]
  
  ## If no non-null tracks, end function.
  if(nrow(tr1) == 0){
    message("Window of available data does not cover time window of interest.")
    message("nrow(tr1) == 0 after removing null tracks")
    return(NULL)
  }
  ## Otherwise, there is at least 1 non-null track:
  # For each sheep, resample to the desired frequency in hours
  # and calculate steps by burst
  
  ## > Resampling ----
  ## Note that "selection is not scale invariant"!!!!
  ## (see Signer et al 2018). For some individuals, sampling is once daily; 
  ## for others it is 2x, 4x, 6x, and 12x daily. This means we will need 
  ## to temporally resample before generating step lengths. Begin by 
  ## removing individuals with infrequent sampling and/or only brief data 
  ## collection windows (frequency and window defined above). 
  
  # Measure sampling rate
  samplingrate = tr1 %>%
    mutate(sr = lapply(trk,
                       summarize_sampling_rate)) %>%
    dplyr::select(ID, sr, Sex) %>%
    unnest(c(sr)) %>%
    mutate(medfreq = lubridate::as.duration(paste0(median," ",unit)))
  samplingrate$durationcalc = (samplingrate$n)*(samplingrate$medfreq)
  retainedInds = samplingrate[samplingrate$mean <= minInt,]$ID
  
  ## Filter to final sheep based on sampling rate
  tr1 = tr1 %>%
    filter(ID %in% retainedInds)
  
  ## If we lose all individuals due to poor sampling frequency, end the function
  if(nrow(tr1) == 0){
    message("Insufficient sampling frequency.")
    message("nrow(tr1) == 0 after removing coarse sampling rate individuals")
    return(NULL)
  }
  # Otherwise, resample to match the desired minimum interval
  # and summarize step length and turning angles
  
  #### Prepare data for random steps ----
  ## Resample track to the pre-defined minimum interval
  tr1resampled <- tr1 %>% 
    mutate(ssf = lapply(trk,
                        function(x){
                          x %>%
                            amt::track_resample(rate = hours(minInt),
                                                tolerance = minutes(15))
                        }))
  ## Measure step lengths and turning angles
  tr1steplength <- tr1resampled %>% 
    mutate(ssf = lapply(ssf,
                        function(x){
                          x %>%
                            steps_by_burst() 
                        }))
  
  ## Finally, return the steplength object.
  return(tr1steplength)
}


## Identify all unique herd unit - year combinations
allHUYrs = dates_df %>%
  mutate(uniqueHuYrID = paste0(HerdUnit,"-",springYr)) %>%
  pull(uniqueHuYrID) %>%
  unique() %>%
  sort()
allHUYrs

## Apply random steps prep by herd unit and year
allSteps = lapply(allHUYrs,
                  randomStepsPrep)

message(paste0("elapsed time: ", (Sys.time() - startTime)))
## *** About 5 minutes

allSteps_list = allSteps[!unlist(lapply(allSteps, is.null))]

allSteps_tibble = do.call("rbind",allSteps_list)
allSteps_tibble$ssf[[1]]

## Combine individual ID and metadata with ssf data
allSteps_metadat = map2(allSteps_tibble$ssf, 
                        allSteps_tibble$ID, 
                        ~.x %>% mutate(ID = .y))
allSteps_metadat = map2(allSteps_metadat, 
                        allSteps_tibble$Year, 
                        ~.x %>% mutate(Year = .y))
allSteps_metadat = map2(allSteps_metadat, 
                        allSteps_tibble$Sex, 
                        ~.x %>% mutate(Sex = .y))
allSteps_metadat = map2(allSteps_metadat, 
                        allSteps_tibble$HU, 
                        ~.x %>% mutate(HU = .y))
allSteps_metadat = map2(allSteps_metadat, 
                        allSteps_tibble$MvmtClass, 
                        ~.x %>% mutate(MvmtClass = .y))
## Convert ssf data with metadata into a single master tibble
ssfMasterTibble = do.call("rbind",
                          allSteps_metadat) %>%
  dplyr::select("ID", "Sex", "HU", "Year", "MvmtClass",
                "burst_", "x1_", "x2_", "y1_", "y2_",
                "sl_", "direction_p", "ta_",
                "t1_", "t2_", "dt_")

#### Generating random steps ----
#### This function takes a population's steps and generates 
#### random steps using the population's distribution of 
#### step lengths and turning angles
makeRandomSteps = function(popStepsObserved){
  #### Generate random steps ----
  ## Generate random steps as background movement possibilities. 
  ## These will serve as a null model for movement decisions.
  ## Movement variability is greatest between midday and night,
  ## so do not include crepuscular differentiation
  popStepsWithRand <- tibble("ssf" = list(popStepsObserved)) %>% 
    mutate(ssf = lapply(ssf,
                        function(x){
                          x %>%
                            time_of_day(include.crepuscule = FALSE, 
                                        where = "both") %>%
                            random_steps(n_control = 30) 
                        }))
  return(popStepsWithRand)
}

#### Execute function ----
myRandomSteps = makeRandomSteps(ssfMasterTibble)

#### Save result ----
saveRDS(myRandomSteps, file = "data/randomSteps.RDS")


#### Epilogue ----
BRRR::skrrrahh(12)
Sys.sleep(5)
stopTime = Sys.time()
message("elapsed time: ")
stopTime-startTime
## Software versions
sessionInfo()
