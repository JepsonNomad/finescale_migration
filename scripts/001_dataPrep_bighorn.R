startTime = Sys.time()
#### Prepping Sierra bighorn movement data for iSSF's
#### Christian John
#### 28 Nov 2023
## Takes about 6.5 min to run

#### Load packages ----
library(plyr)
library(lubridate)
library(ggplot2)
library(RStoolbox)
library(tidyverse)
library(sf)
library(rgdal)
library(terra)
library(amt)
library(purrr)
library(migrateR)
library(adehabitatLT)
source("scripts/000_migrateRoverrides.R")
source("scripts/000_samplingParameters.R")


#### Import data ----
## > Herd Units ----
HU_raw = read_sf("../SNBS_HerdUnits/Cathedral/snbs_herdunits_Cathedral_BlkDvd.shp")
HU_buf = read_sf("../SNBS_HerdUnits/Buffered/SNBS_HUs_5kmUnion.shp") %>%
  st_transform(crs(HU_raw))


## > Movement data ----
# Import movement data 
sheep_imported = read_csv("../data_ts/SNBS_Locs/SNBS_Herd_Sex_20231115_KAClean.csv")

spatialOutliers = sheep_imported %>%
  filter(Longitude > -80)
## Remove the faulty Maryland data point
## Fix missing herd units and sex for S182, S508, and S213
sheep_selected = sheep_imported %>% 
  filter(Longitude < -80) %>%
  mutate(Herd_Rels = case_when(Animal_ID == "S182" ~ "Wh",
                               Animal_ID == "S508" ~ "Wh",
                               Animal_ID == "S213" ~ "Bx",
                               TRUE ~ Herd_Rels)) %>%
  mutate(Sheep_Sex = case_when(Animal_ID == "S182" ~ "F",
                               Animal_ID == "S508" ~ "M",
                               Animal_ID == "S213" ~ "F",
                               TRUE ~ Sheep_Sex))
sheep_sf = st_as_sf(sheep_selected,
                    coords = c("Longitude", "Latitude"),
                    crs = 4326) %>%
  mutate(dateTime = lubridate::mdy_hms(DateTimePST))

s192 = sheep_sf %>%
  filter(Animal_ID == "S192") %>%
  st_transform(crs(HU_raw))
s192

ggplot() +
  geom_sf(data = HU_buf, fill = "grey40") +
  geom_sf(data = HU_raw, fill = "grey80", col = "orange2") +
  geom_sf(data = s192, aes(col = as.Date(dateTime))) +
  scale_color_viridis_c(option="E", trans = "date")

## > DEM ----
## Fine-scale DEM
Terrain = terra::rast("../USGS_3DEP/data/DEM_UTM.tif")

# Check that movement data are projected the same as habitat data
if(st_crs(sheep_sf) != crs(Terrain)){
  sheep_sf = sheep_sf %>%
    st_transform(st_crs(Terrain))
}


#### Data wrangling ----
## > DEM ----
# Crop DEM to extent of sheep locations
Terrain = crop(Terrain, vect(sheep_sf))

## > Spatiotemporal bighorn data ----
sheep_raw = sheep_sf %>%
  mutate(UTME = st_coordinates(.)[,1],
         UTMN = st_coordinates(.)[,2]) %>%
  st_drop_geometry()

firstLoc = min(sheep_raw$dateTime)
# [1] "2002-03-18 02:15:36 UTC"
lastLoc = max(sheep_raw$dateTime)
# [1] "2023-10-30 00:00:44 UTC"
## Wrangle names
sheep_raw$x = sheep_raw$UTME
sheep_raw$y = sheep_raw$UTMN
sheep_raw$t = sheep_raw$dateTime

## > Bighorn individuals ----
## Format for multiple individuals (list-style)
# Generate list of data for multiple individual analysis
dat_all = sheep_raw %>%
  nest(-Animal_ID) %>%
  rename(ID = Animal_ID)

## Add individual-level parameters to list
# Extract sex from raw data and assign to list
Sex = sheep_raw %>%
  group_by(Animal_ID) %>%
  summarize(IDsex = unique(Sheep_Sex))
dat_all$Sex = Sex$IDsex[match(dat_all$ID, Sex$Animal_ID)]
# Extract herd unit from raw data and assign to list
HU = sheep_raw %>%
  group_by(Animal_ID) %>%
  summarize(IDHU = unique(Herd_Rels))
dat_all$HU = HU$IDHU[match(dat_all$ID, HU$Animal_ID)]

# Reorganize because it's prettier to look at
dat_all = dat_all %>%
  dplyr::select(ID, HU, Sex, data)

#### amt data prep ----
## > Getting into the amt-verse ----
# Convert to a track object
tr1 = dat_all %>%
  mutate(trk = lapply(data,
                      function(d){
                        amt::make_track(d,x,y,t,
                                        crs = st_crs(sheep_sf))
                      }))
tail(tr1)


## > Filter by sampling parameters ----
# Retain only individuals meeting minimum sampling requirements (defined at top)
# Note that there are 3 cases with fewer than 10 total observations.
# These obviously aren't useful in our analyses of mvmt @ seasonal level so get 
# rid of them now.
# Measure sampling rate
samplingrate = tr1 %>%
  filter(lapply(trk, function(x){nrow(x)}) > 10) %>%
  mutate(sr = lapply(trk,
                     summarize_sampling_rate)) %>%
  dplyr::select(ID, sr, Sex) %>%
  unnest(c(sr)) %>%
  mutate(medfreq = lubridate::as.duration(paste0(median," ",unit)))
samplingrate$durationcalc = (samplingrate$n)*(samplingrate$medfreq)

# Identify individuals to remove, based on:
# sampling frequency...
REMOVEfreq = samplingrate %>%
  filter(medfreq > lubridate::as.duration(minFreq)) %>%
  dplyr::pull(ID) %>%
  as.character()
# ...and sampling duration
REMOVEdura = samplingrate %>%
  filter(durationcalc < as.duration(minDura)) %>%
  dplyr::pull(ID) %>%
  as.character()
REMOVEINDS = unique(c(REMOVEfreq, REMOVEdura))

# Remove individuals that did not meet sampling requirements
tr1 = tr1 %>%
  filter(!(as.character(ID) %in% REMOVEINDS))
message(paste0("Retaining ", nrow(tr1)," individuals"))
## Retaining 331 individuals
retainedInds = tr1$ID

# Create spatial object with only retained individual bighorn data
retained_sf = sheep_sf[as.character(sheep_sf$Animal_ID) %in% retainedInds,]

# Extract elevation for bighorn locations
retained_sf$elev = extract(Terrain, vect(retained_sf))[,2]


#### Classifying migration ----
## > Define functions ----
## A function to bin one sheep-year of movement into migration classes
# Classify individuals' movement patterns using
# a minimum residency period (rho) and distance between 
# ranges (delta) based on parameters defined at top.
mvmtClassifier = function(burstID, 
                          dat, 
                          crsTarget, 
                          min.dura = minDura,
                          min.rho = min.r, 
                          min.delta = min.d){
  subind = dat[dat$indBurst == burstID,]
  indHU = subind$Herd_Rels[1]
  # Reformat starting date to account for leap years
  fallYr = str_split(burstID,pattern="-")[[1]][1]
  migRstdt = as.Date(paste0(fallYr,"-",startDateMigR),
                     format = "%Y-%j")
  migR_monDay = strftime(migRstdt,
                         format = "%m-%d")
  
  minDays = as.numeric(str_split(min.dura,pattern=" ")[[1]][1])*7
  ## Only analyze if the spring data make it through February at least
  ## and for individuals with at least 50 unique observation days 
  ## Note min.dura*7 = the range of dates because minDura was defined in weeks
  if(length(unique(subind$Date)) >= minDays){
    firstDat = as.Date(min(subind$dateTime))
    subind_sf = st_as_sf(subind,
                         coords = c("UTME", "UTMN"),
                         crs = crsTarget)
    ind.ltraj <- as.ltraj(st_coordinates(subind_sf),
                          date=as.POSIXct(subind$dateTime), 
                          id=subind$Animal_ID[1],
                          infolocs = data.frame(elev = 
                                                  subind$elev))
    ## Name of infolocs MUST BE "elev" and for 
    ## whatever reason adehabitatLT makes it elevation
    for(i in 1:length(ind.ltraj)){
      names(infolocs(ind.ltraj)[[i]]) <- "elev"
    }
    
    ## Apply migrateR workflow to the ltraj object
    subind.mvmt = mvmtClass(ind.ltraj, 
                            warnOnly = TRUE,
                            fam = "elev", 
                            p.est = pEst(s.d = -1500,
                                         u.d = 0))
    subind.mvmt.ref = refine(subind.mvmt, 
                             p.est = pEst(s.d = -1500))
    subind.topmvmt = topmvmt(subind.mvmt.ref, 
                             mdelta = min.delta,
                             mrho = min.rho)
    # Format the data to suit data tidying
    subind.list = mvmt2df(subind.topmvmt)
    # Get the type of model that was selected
    topModel = names(subind.list)
    if(topModel == "migrant"){
      t2 = theta2(subind.mvmt.ref, topModel)
      subind.list$migrant["theta2"] <- unlist(t2["theta2"])
    }
    subind.df = as.data.frame(subind.list[[topModel]])
    subind.df$burstYr = burstID
    subind.df$mvmtClass = topModel
    subind.df$indYrStartDate = firstDat
    subind.df$ID = as.character(unique(subind$Animal_ID))[1]
    subind.df$Sex = unique(subind$Sheep_Sex)
    subind.df$HerdUnit = indHU
    return(subind.df)
  }else{
    return(data.frame(gamma = NA,
                      theta = NA,
                      theta2 = NA,
                      phi = NA,
                      phi2 = NA,
                      delta = NA,
                      burstYr = burstID,
                      mvmtClass = NA,
                      ID = unique(subind$Animal_ID),
                      Sex = unique(subind$Sheep_Sex),
                      HerdUnit = indHU))   
  }
}

## A function to route an individual sheep-year into the migrateR workflow
# This function takes a single individual and a pre-defined annual breakpoint
# Divides that individual's GPS data into years bounded by the breakpoint
# And applies the mvmtClassifier to each sheep-year for that individual
SheepAnalysis = function(indID,
                         yrBreakPt = startDateMigR){
  ## First, set a few custom breakpoints to deal with animals that undertook 
  ## unusually timed movements
  if(indID == "S105"){
    ## S105 model settles on weirdly 'late migrations' because of late 
    ## season movements from mid elevs to high
    yrBreakPt <- 244
  }
  if(indID == "S251"){
    ## For a cleaner symmetry around seasonal movements
    yrBreakPt <- 244
  }
  if(indID == "S295"){
    ## S295 weird early winter movement. revise breakpt
    yrBreakPt <- 273
  }
  if(indID == "S495"){
    ## For a cleaner symmetry around movements
    yrBreakPt <- 273
  }
  if(indID == "S526"){
    ## For a cleaner symmetry around movements
    yrBreakPt <- 273
  }
  if(indID == "S575"){
    ## Work around a quick winter foray
    yrBreakPt <- 244
  }
  
  # Create one individual's worth of individual annual GPS location data
  # (i.e. one sheep-year) and create some useful date/time columns
  mySheep = retained_sf[retained_sf$Animal_ID==indID,] %>%
    mutate(Date = as.Date(dateTime),
           DOY = yday(dateTime))
  
  # Take out duplicate observations
  if(sum(duplicated(mySheep$Date))>0){
    mySheep<-mySheep[which(!duplicated(mySheep$Date)),] 
  }
  # Take out NA dates
  if(sum(is.na(mySheep$Date))>0){
    mySheep = mySheep[(-which(is.na(mySheep$Date))),] 
  }
  # Convert to a plain old data.frame()
  ind.df <- mySheep %>%
    mutate(UTME = st_coordinates(.)[,1],
           UTMN = st_coordinates(.)[,2]) %>% 
    st_drop_geometry()
  ind.df = ind.df %>%
    mutate(BreakSlot = ifelse(DOY < yrBreakPt,
                              "Pre",
                              "Post"),
           yrBreak = ifelse(BreakSlot == "Post",
                            as.numeric(Year) + 1,
                            as.numeric(Year)),
           indBurst = paste0(yrBreak-1,"-",yrBreak))
  
  ## Apply the movement classifier to all individual sheep-years
  ## for this individual
  ind.bursts = unique(ind.df$indBurst)
  ind.subdfs = lapply(ind.bursts,
                      mvmtClassifier,
                      crsTarget = st_crs(mySheep),
                      dat = ind.df)
  ## Find date ranges for each individual sheep-year
  ind.subdfs = lapply(ind.subdfs,
                      function(x){
                        b = x$burstYr
                        x$indYrFirstObs = min(
                          as.Date(
                            ind.df$Date[ind.df$indBurst==b]))
                        x$indYrLastObs = max(
                          as.Date(
                            ind.df$Date[ind.df$indBurst==b]))
                        return(x)
                      })
  ## Combine
  ind.mvmtparams = do.call("rbind.fill",
                           ind.subdfs) 
  ## If no years' worth of data were useful, ind.mvmtparams will have very few columns
  ## If only resident classes were detected, no theta col will be present
  ## Only when dispersers/migrants are present will theta values be available
  if(sum(is.na(ind.mvmtparams$mvmtClass)) == nrow(ind.mvmtparams)){
    ind.mvmtparams = ind.mvmtparams %>%
      separate(burstYr, into = c("fallYr","springYr"),
               remove=FALSE)
  }else if(sum(ind.mvmtparams$mvmtClass=="resident",
               na.rm = TRUE) == nrow(ind.mvmtparams)){
    ind.mvmtparams = ind.mvmtparams %>%
      separate(burstYr, into = c("fallYr","springYr"),
               remove=FALSE) %>%
      mutate(indYrStartDOY = strftime(indYrStartDate, "%j"),
             indYrEndYr = ifelse(indYrStartDOY < yrBreakPt,
                                 as.numeric(strftime(indYrStartDate, 
                                                     "%Y")),
                                 as.numeric(strftime(indYrStartDate + 
                                                       years(1), 
                                                     "%Y"))),
             indYrEndDate = as.Date(paste0(indYrEndYr,"-",yrBreakPt),
                                    format = "%Y-%j"),
             startElev = gamma)
  }else if(sum(ind.mvmtparams$mvmtClass != "migrant", na.rm=T) == nrow(ind.mvmtparams)){
    ind.mvmtparams = ind.mvmtparams %>%
      separate(burstYr, into = c("fallYr","springYr"),
               remove=FALSE) %>%
      mutate(indYrStartDOY = strftime(indYrStartDate, "%j"),
             indYrEndYr = ifelse(indYrStartDOY < yrBreakPt,
                                 as.numeric(strftime(indYrStartDate, 
                                                     "%Y")),
                                 as.numeric(strftime(indYrStartDate + 
                                                       years(1), 
                                                     "%Y"))),
             indYrEndDate = as.Date(paste0(indYrEndYr,"-",yrBreakPt),
                                    format = "%Y-%j"),
             migTiming1 = indYrStartDate + days(round(theta)),
             startElev = gamma,
             endElev = gamma+delta)
  }else{
    ind.mvmtparams = ind.mvmtparams %>%
      separate(burstYr, into = c("fallYr","springYr"),
               remove=FALSE) %>%
      mutate(indYrStartDOY = strftime(indYrStartDate, "%j"),
             indYrEndYr = ifelse(indYrStartDOY < yrBreakPt,
                                 as.numeric(strftime(indYrStartDate, 
                                                     "%Y")),
                                 as.numeric(strftime(indYrStartDate + 
                                                       years(1), 
                                                     "%Y"))),
             indYrEndDate = as.Date(paste0(indYrEndYr,"-",yrBreakPt),
                                    format = "%Y-%j"),
             migTiming1 = indYrStartDate + days(round(theta)),
             migTiming2 = indYrStartDate + days(round(theta2)),
             startElev = gamma,
             endElev = gamma+delta)
  }
  return(ind.mvmtparams)
}

## > Apply functions ----
## Analyze each individual-year and recombine the data
migrSummary = lapply(retainedInds,
                     SheepAnalysis)
HU_migrations = do.call("rbind.fill",
                        migrSummary)
migrDatesDF = HU_migrations
# Sometimes the derivation of theta2 generates out-of-sample metrics; 
# in that case, remove them and treat movement as disperser
migrDatesDF$mvmtClass[migrDatesDF$migTiming2 > migrDatesDF$indYrLastObs] <- "disperser"
migrDatesDF$migTiming2[migrDatesDF$migTiming2 > migrDatesDF$indYrLastObs] <- NA


#### Save results ----
write_csv(migrDatesDF,
          file = paste0("data/migrationPatterns_allInds.csv"))


#### Epilogue ----
BRRR::skrrrahh(12)
Sys.sleep(0.5)
print(sessionInfo())
stopTime = Sys.time()
print(stopTime-startTime)

fileConn<-file("sessInfo/001_dataPrep.txt")
writeLines(capture.output(print(stopTime-startTime)),fileConn)
writeLines("\n\n\n",fileConn)
writeLines(capture.output(sessionInfo()),fileConn)
close(fileConn)
