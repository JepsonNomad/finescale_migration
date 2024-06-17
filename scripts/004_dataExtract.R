startTime = Sys.time()

## Use Sierra bighorn steps and random steps to extract
## habitat covariates for use in step selection functions
## ~11 minutes collectively

#### Import packages ----
library(lubridate)
library(RStoolbox)
library(tidyverse)
library(sf)
library(rgdal)
library(terra)
library(amt)
library(purrr)
source("scripts/000_extract_covariates.R")


#### Define parameters ----
outBase = "data/extract"
if(!dir.exists(outBase)){
  dir.create(outBase)
}

set.seed(1927) # Invention of Pez

#### Load data ----
## Movement data
# Choose crs
sheepCRS = st_crs(32611)

# Get selected and random steps from 003_randomSteps.R
myRandomSteps = readRDS("data/randomSteps.rds")


#### Write function ----
#### IMPERATIVE: Ensure that crs is consistent across data!!!
#### Through random steps generation, we have lost the crs
#### information for sheep location data.
#### This function takes a population's steps and random steps,
#### extracts all covariates (static and temporally variant)
#### and saves data to a csv.
## mySteps should be a tibble with a column called ssf, which contains a list of ssf tibble data
randStepsExtract = function(mySteps = myRandomSteps,
                            stepsCRS = sheepCRS,
                            myYear = 2005,
                            outFolder = outBase){
  #### Identify dates of interest ----
  myYearTibble = mySteps$ssf[[1]] %>%
    filter(Year == myYear)
  mySSF = tibble("ssf" = list(myYearTibble))
  selectedCoords = myYearTibble %>%
    filter(case_)
  randomCoords = myYearTibble %>%
    filter(!case_)
  myDateRangeTibble = myYearTibble %>%
    dplyr::select(t1_,t2_) %>%
    summarize(minDate = min(t1_),
              maxDate = max(t2_))
  thisYrMinDate = as.Date(pull(myDateRangeTibble[1,1]))
  thisYrMaxDate = as.Date(pull(myDateRangeTibble[1,2]))
  
  #### BEGIN Data loading component ----
  ## .... Fractional snow cover ----
  SSN_FSC_FPlist1 = list.files(path = paste0(
    "../SSN/data/",
    as.numeric(myYear)-1),
    pattern = "*.tif",
    full.names = TRUE)
  SSN_FSC_FPlist2 = list.files(path = paste0(
    "../SSN/data/",
    myYear),
    pattern = "*.tif",
    full.names = TRUE)
  SSN_FSC_FPlist = c(SSN_FSC_FPlist1,SSN_FSC_FPlist2)
  # SSN_FSC_FPlist
  SSN_FSC_FPlist_dates = sapply(str_split(SSN_FSC_FPlist,
                                          "[.]"),
                                "[[",
                                3) %>%
    as.Date(format="%Y%m%d")
  FSC_datesOfInterestIndices = which(SSN_FSC_FPlist_dates >= thisYrMinDate &
                                       SSN_FSC_FPlist_dates <= thisYrMaxDate)
  selectedFSCDates = SSN_FSC_FPlist_dates[FSC_datesOfInterestIndices]
  FSC = terra::rast(SSN_FSC_FPlist[FSC_datesOfInterestIndices])
  terra::time(FSC) <- selectedFSCDates
  
  ## .... Snow phenology ----
  SSN_snowphenoFP = paste0(
    "../SSN/snowPhenologyMosaics/",
    myYear,
    ".tif")
  snowPhen = terra::rast(SSN_snowphenoFP)
  names(snowPhen) <- c("fallMidpoint","meltMidpoint","fallRate","meltRate")
  
  ## .... Distance to snow ----
  ## List layer files of interest
  SSN_snowDistanceFPs1 = list.files(paste0(
    "../SSN/data/",
    as.numeric(myYear)-1,
    "/snowDistance"),
    pattern = "*.tif",
    full.names = TRUE)
  SSN_snowDistanceFPs2 = list.files(paste0(
    "../SSN/data/",
    myYear,
    "/snowDistance"),
    pattern = "*.tif",
    full.names = TRUE)
  SSN_snowDistanceFPs = c(SSN_snowDistanceFPs1,SSN_snowDistanceFPs2)
  ## Find dates of snow distance layers
  SSN_snowDistanceFP_dates = sapply(
    str_split(SSN_snowDistanceFPs,
              "[.]"),
    "[[",
    3) %>%
    as.Date(format="%Y%m%d")
  ## Subset snow layers based on dates of interest
  snowDist_datesOfInterestIndices = which(
    SSN_snowDistanceFP_dates >= thisYrMinDate &
      SSN_snowDistanceFP_dates <= thisYrMaxDate)
  selectedSnowDistDates =
    SSN_snowDistanceFP_dates[snowDist_datesOfInterestIndices]
  snowDist = terra::rast(
    SSN_snowDistanceFPs[snowDist_datesOfInterestIndices])
  
  # Rename snow distance layers
  rawDistNames_reform = as.character(format(selectedSnowDistDates,
                                            format = "%Y%m%d"))
  snowDist_layerNames = rawDistNames_reform
  # Assign timestamp to MODIS layers
  terra::time(snowDist) <- selectedSnowDistDates
  # snowDist <- raster::setZ(snowDist, 
  #                          selectedSnowDates)
  
  ## .... NDVI ----
  ## Plant phenology
  MOD_phenoFP = paste0(
    "../MOD13Q1/data/QA1/phenology_",
    myYear,
    ".tif")
  ndviPhen = terra::rast(MOD_phenoFP)
  names(ndviPhen) <- c("greenupMidpoint","greenupRate",
                       "senescMidpoint","senescRate")
  MOD_scaleFP = paste0(
    "../MOD13Q1/data/QA1/scaleLims_",
    myYear,
    ".tif")
  ndviScal = terra::rast(MOD_scaleFP)
  names(ndviScal) <- c("wNDVI","mNDVI")
  
  ## .... Daymet ----
  ## Identify the layers of interest
  DAYMET_layerDates1 = read.csv(
    "../DAYMET/data/compiled/layerDates.csv") %>%
    filter(Year == as.numeric(myYear)-1) %>%
    mutate(Date = as.Date(paste0(Year,"-",Month,"-",Day)),
           dateChar = as.character(Date)) %>%
    filter(Date >= thisYrMinDate,
           Date <= thisYrMaxDate)
  DAYMET_layerDates2 = read.csv(
    "../DAYMET/data/compiled/layerDates.csv") %>%
    filter(Year == myYear) %>%
    mutate(Date = as.Date(paste0(Year,"-",Month,"-",Day)),
           dateChar = as.character(Date)) %>%
    filter(Date >= thisYrMinDate,
           Date <= thisYrMaxDate)
  DAYMET_layerDates = rbind(DAYMET_layerDates1,
                            DAYMET_layerDates2)
  ## Import the appropriate layers
  highTemps = terra::rast(
    "../DAYMET/data/compiled/tmax.tif",
    lyrs = DAYMET_layerDates$layerIndex)
  precip = terra::rast(
    "../DAYMET/data/compiled/precip.tif",
    lyrs = DAYMET_layerDates$layerIndex)
  
  terra::time(highTemps) <- as.Date(DAYMET_layerDates$dateChar)
  terra::time(precip) <- as.Date(DAYMET_layerDates$dateChar)
  
  ## .... DEM ----
  Terrain = terra::rast(list.files("../USGS_3DEP/data",
                                   pattern = "*.tif",
                                   full.names = TRUE))
  names(Terrain) = c("ASP_UTM","DEM_UTM",
                     "ESC_Distance_UTM",
                     "ESC_UTM","SLP_UTM")
  
  ## .... Landcover ----
  NLCD = terra::rast("../NLCD/Counties/NLCD_SNBS_counties_UTM.tif")
  names(NLCD) <- "NLCD"
  #### END data loading component ----
  message("..Datasets loaded.")
  
  
  ## > CRS check ----
  # Check that movement data are projected the same as habitat data
  if(stepsCRS != st_crs(FSC)){
    message("WARNING: SNOW COVER DATA NOT IN SAME PROJECTION AS SHEEP DATA!")
  }
  if(stepsCRS != st_crs(snowPhen)){
    message("WARNING: SNOW PHENOLOGY NOT IN SAME PROJECTION AS SHEEP DATA!")
  }
  if(stepsCRS != st_crs(snowDist)){
    message("WARNING: SNOW DISTANCE NOT IN SAME PROJECTION AS SHEEP DATA!")
  }
  if(stepsCRS != st_crs(ndviPhen)){
    message("WARNING: NDVI PHENOLOGY DATA NOT IN SAME PROJECTION AS SHEEP DATA!")
  }
  if(stepsCRS != st_crs(highTemps)){
    message("WARNING: TEMPERATURE NOT IN SAME PROJECTION AS SHEEP DATA!")
  }
  if(stepsCRS != st_crs(precip)){
    message("WARNING: PRECIPITATION NOT IN SAME PROJECTION AS SHEEP DATA!")
  }
  if(stepsCRS != st_crs(Terrain)){
    message("WARNING: TERRAIN NOT IN SAME PROJECTION AS SHEEP DATA!")
  }
  if(stepsCRS != st_crs(NLCD)){
    message("WARNING: LAND COVER NOT IN SAME PROJECTION AS SHEEP DATA!")
  }
  
  #### Crop covariate layers ----
  ## Trim covariates to ROI
  xmn = min(c(myYearTibble$x1_,myYearTibble$x2_)-10)
  xmx = max(c(myYearTibble$x1_,myYearTibble$x2_)+10)
  ymn = min(c(myYearTibble$y1_,myYearTibble$y2_)-10)
  ymx = max(c(myYearTibble$y1_,myYearTibble$y2_)+10)
  
  ## Create an appropriate extent
  r <- rast()
  ext(r) <- c(xmn,xmx,ymn,ymx)
  e <- ext(r)
  
  ## Crop covariate layers to extent
  # FSC = terra::crop(FSC,e)
  # snowDist = terra::crop(snowDist,e)
  # snowPhen = terra::crop(snowPhen,e)
  # ndviPhen = terra::crop(ndviPhen,e)
  # ndviScal = terra::crop(ndviScal,e)
  # highTemps = terra::crop(highTemps,e)
  # precip = terra::crop(precip,e)
  # Terrain = terra::crop(Terrain,e)
  # NLCD = terra::crop(NLCD,e)
  
  message("..Covariate layers cropped.")
  
  #### Prep movement data ----
  #### Data extraction ----
  ## Dynamic layers (variable through time)
  ## Fractional snow cover
  message("...extracting snow cover")
  tr1covsFSC <- mySSF %>%
    mutate(ssf = lapply(ssf,
                        function(x){
                          x %>%
                            extract_covariates_var_time(
                              FSC,
                              max_time = days(1), 
                              where = "both") %>%
                            mutate(FSC_start = time_var_covar_start,
                                   FSC_end = time_var_covar_end)
                        }))
  ## Distance to snow
  message("...extracting snow distance")
  tr1covsSnowDist <- mySSF %>%
    mutate(ssf = lapply(ssf,
                        function(x){
                          x %>%
                            extract_covariates_var_time(
                              snowDist,
                              max_time = days(1), 
                              where = "both") %>%
                            mutate(
                              snowDist_start = time_var_covar_start,
                              snowDist_end = time_var_covar_end)
                        }))
  ## High temperature
  message("...extracting temperature")
  tr1covsHighTemp <- mySSF %>%
    mutate(ssf = lapply(ssf,
                        function(x){
                          x %>%
                            extract_covariates_var_time(
                              highTemps,
                              max_time = days(1), 
                              where = "both") %>%
                            mutate(
                              highTemp_start = time_var_covar_start,
                              highTemp_end = time_var_covar_end)
                        }))
  ## Precipitation
  message("...extracting precipitation")
  tr1covsPrecip <- mySSF %>%
    mutate(ssf = lapply(ssf,
                        function(x){
                          x %>%
                            extract_covariates_var_time(
                              precip,
                              max_time = days(1), 
                              where = "both") %>%
                            mutate(precip_start = time_var_covar_start,
                                   precip_end = time_var_covar_end)
                        }))
  
  ## Static layers (static through time)
  ## Note that the snow phenology layer will be used dynamically
  ## when we calculate deltaMelt
  ## Snow phenology
  message("...extracting snow phenology")
  tr1covsSnowPhen <- mySSF %>%
    mutate(ssf = lapply(ssf,
                        function(x){
                          x %>%
                            amt::extract_covariates(
                              snowPhen, 
                              where = "both") 
                        }))
  ## NDVI phenology
  message("...extracting ndvi phenology")
  tr1covsNDVIPhen <- mySSF %>%
    mutate(ssf = lapply(ssf,
                        function(x){
                          x %>%
                            amt::extract_covariates(
                              ndviPhen, 
                              where = "both") 
                        }))
  ## NDVI scale limits
  message("...extracting ndvi scale limits")
  tr1covsNDVIScal <- mySSF %>%
    mutate(ssf = lapply(ssf,
                        function(x){
                          x %>%
                            amt::extract_covariates(
                              ndviScal, 
                              where = "both") 
                        }))
  
  ## Terrain
  message("...extracting terrain")
  tr1covsTerrain <- mySSF %>%
    mutate(ssf = lapply(ssf,
                        function(x){
                          x %>%
                            extract_covariates(
                              Terrain, 
                              where = "both") 
                        }))
  # print(Sys.time())
  ## Land cover
  message("...extracting land cover")
  tr1covsNLCD <- mySSF %>%
    mutate(ssf = lapply(ssf,
                        function(x){
                          x %>%
                            extract_covariates(
                              NLCD, 
                              where = "both") 
                        }))
  # print(Sys.time())

  message("...compiling extracted covariates")
  #### Combine results ----
  ## Recompile all extracted covariates into one ssf dataframe
  tr1covsCompiled = tr1covsSnowPhen
  ## NDVI phenology (MODIS)
  tr1covsCompiled$ssf[[1]]$greenupMidpoint_start =
    tr1covsNDVIPhen$ssf[[1]]$greenupMidpoint_start
  tr1covsCompiled$ssf[[1]]$greenupMidpoint_end =
    tr1covsNDVIPhen$ssf[[1]]$greenupMidpoint_end
  
  tr1covsCompiled$ssf[[1]]$greenupRate_start =
    tr1covsNDVIPhen$ssf[[1]]$greenupRate_start
  tr1covsCompiled$ssf[[1]]$greenupRate_end =
    tr1covsNDVIPhen$ssf[[1]]$greenupRate_end
  
  tr1covsCompiled$ssf[[1]]$senescMidpoint_start =
    tr1covsNDVIPhen$ssf[[1]]$senescMidpoint_start
  tr1covsCompiled$ssf[[1]]$senescMidpoint_end =
    tr1covsNDVIPhen$ssf[[1]]$senescMidpoint_end
  
  tr1covsCompiled$ssf[[1]]$senescRate_start =
    tr1covsNDVIPhen$ssf[[1]]$senescRate_start
  tr1covsCompiled$ssf[[1]]$senescRate_end =
    tr1covsNDVIPhen$ssf[[1]]$senescRate_end
  
  ## NDVI double logistic scaling values
  tr1covsCompiled$ssf[[1]]$mNDVI_start =
    tr1covsNDVIScal$ssf[[1]]$mNDVI_start
  tr1covsCompiled$ssf[[1]]$mNDVI_end =
    tr1covsNDVIScal$ssf[[1]]$mNDVI_end
  
  tr1covsCompiled$ssf[[1]]$wNDVI_start =
    tr1covsNDVIScal$ssf[[1]]$wNDVI_start
  tr1covsCompiled$ssf[[1]]$wNDVI_end =
    tr1covsNDVIScal$ssf[[1]]$wNDVI_end
  
  ## Fractional snow cover (SSN)
  tr1covsCompiled$ssf[[1]]$FSC_start =
    tr1covsFSC$ssf[[1]]$FSC_start
  tr1covsCompiled$ssf[[1]]$FSC_end = 
    tr1covsFSC$ssf[[1]]$FSC_end
  ## Snow melt phenology (SSN)
  tr1covsCompiled$ssf[[1]]$snowDist_start =
    tr1covsSnowDist$ssf[[1]]$snowDist_start
  tr1covsCompiled$ssf[[1]]$snowDist_end = 
    tr1covsSnowDist$ssf[[1]]$snowDist_end
  ## High tempeature (DAYMET)
  tr1covsCompiled$ssf[[1]]$highTemp_start =
    tr1covsHighTemp$ssf[[1]]$highTemp_start
  tr1covsCompiled$ssf[[1]]$highTemp_end =
    tr1covsHighTemp$ssf[[1]]$highTemp_end
  ## Precipitation (DAYMET)
  tr1covsCompiled$ssf[[1]]$Precip_start = 
    tr1covsPrecip$ssf[[1]]$precip_start
  tr1covsCompiled$ssf[[1]]$Precip_end = 
    tr1covsPrecip$ssf[[1]]$precip_end
  ## Terrain - elevation
  tr1covsCompiled$ssf[[1]]$DEM_start = 
    tr1covsTerrain$ssf[[1]]$DEM_UTM_start
  tr1covsCompiled$ssf[[1]]$DEM_end = 
    tr1covsTerrain$ssf[[1]]$DEM_UTM_end
  ## Terrain - slope
  tr1covsCompiled$ssf[[1]]$SLP_start = 
    tr1covsTerrain$ssf[[1]]$SLP_UTM_start
  tr1covsCompiled$ssf[[1]]$SLP_end = 
    tr1covsTerrain$ssf[[1]]$SLP_UTM_end
  ## Terrain - aspect
  tr1covsCompiled$ssf[[1]]$ASP_start = 
    tr1covsTerrain$ssf[[1]]$ASP_UTM_start
  tr1covsCompiled$ssf[[1]]$ASP_end = 
    tr1covsTerrain$ssf[[1]]$ASP_UTM_end
  ## Terrain - dist from escape terrain
  tr1covsCompiled$ssf[[1]]$ESC_Dist_start =
    tr1covsTerrain$ssf[[1]]$ESC_Distance_UTM_start
  tr1covsCompiled$ssf[[1]]$ESC_Dist_end =
    tr1covsTerrain$ssf[[1]]$ESC_Distance_UTM_end
  ## NLCD
  tr1covsCompiled$ssf[[1]]$NLCD_start =
    tr1covsNLCD$ssf[[1]]$NLCD_start
  tr1covsCompiled$ssf[[1]]$NLCD_end =
    tr1covsNLCD$ssf[[1]]$NLCD_end
  
  tr1Covariates = tr1covsCompiled$ssf[[1]] %>%
    mutate(stepID = paste0(ID,"_",Year,"_",step_id_))
  
  message("..Covariate data compiled.")
  
  #### Save and return result ----
  ## Save output
  write_csv(tr1Covariates,
            file = paste0(outFolder,"/cov_extract_",myYear,".csv"))
  ## Return values
  return(tr1Covariates)
}



#### Execute function ----
allYrs = myRandomSteps$ssf[[1]] %>%
  filter(Year >= 2003,
         Year <= 2022) %>%
  pull(Year) %>%
  unique() %>%
  sort()
i=allYrs[1]
for(i in allYrs){
  message(i)
  strtTime = Sys.time()
  message(".Start time: ")
  message(paste0("..",strtTime))
  randStepsExtract(myYear = i)
  stopTime = Sys.time()
  message(".Elapsed time: ")
  message(paste0("..",(stopTime-strtTime)))
}


#### Epilogue ----
stopTime = Sys.time()
message("Total elapsed time: ")
stopTime-startTime
BRRR::skrrrahh(12)
Sys.sleep(5)
## Software versions
sessionInfo()

