startTime = Sys.time()

## A script to compile extracted habitat covariates for
## step selection modeling
## Takes about 3.5 minutes

#### Import packages ----
library(tidyverse)
library(amt)

## The starting month-day for snow phenology calculations
startMonthDay = "08-15"

#### Define functions ----
## Double logistic NDVI function
ndviCurve = function(DOY, xmidS, xmidA, scalS, scalA){
  NDVI = 
    ((1/(1+exp((xmidS-DOY)/scalS))) -
       (1/(1+exp((xmidA-DOY)/scalA))))
  return(NDVI)
}
## A function to check for high-variance individuals
highVarInds = function(mySSF){
  # Quick check for ssf fit issues
  # Flag high-variance individuals
  summarizeSSF = list()
  for(i in 1:nrow(mySSF)){
    indSSF = mySSF[i,]
    indMvmtRate = indSSF$uphillRate
    ind = as.character(indSSF$ID)
    indMod = indSSF$model[[1]]
    indSum = summary(indMod)
    indCoe = indSum$coefficients
    indParam = rownames(indCoe)
    indSlope = indCoe[,1]
    indSE = indCoe[,3]
    summarizeSSF$ID[[i]] = ind
    summarizeSSF$upRate[[i]] = indMvmtRate
    summarizeSSF$Slope[[i]] = data.frame("param" = indParam,
                                         "slope" = indSlope)
    summarizeSSF$SE[[i]] = data.frame("param" = indParam,
                                      "SE" = indSE)
  }
  for(i in 1:length(summarizeSSF$ID)){
    y = summarizeSSF$SE[i][[1]]
    ID = summarizeSSF$ID[i]
    upRate = summarizeSSF$upRate[i]
    highSEs = y[y$SE > 2,]
    if(nrow(highSEs) > 0){
      print(paste0("***** INDIVIDUAL *****: ",
                   as.character(ID)))
      print(paste0("High SE detected for: ",
                   as.character(highSEs$param)))
      print(paste0("SE_",
                   as.character(highSEs$param),
                   " = ",
                   as.character(highSEs$SE)))
    }
  }
}
## A plotting function for population-level coefficients
plot_popCoeffs = function(mySSF,
                          drawplot=TRUE){
  # First generate overall summaries for coefficient estimates
  d1 <- mySSF %>% 
    mutate(coef = map(model, ~ broom::tidy(.x$model))) %>%
    select(ID, Sex, uphillRate, mvmtClass, coef) %>% 
    unnest(cols = c(coef)) %>%
    mutate(ID = factor(ID)) %>% 
    group_by(term) %>%
    summarize(
      mean = mean(estimate),
      ymin = mean - 1.96 * sd(estimate),
      ymax = mean + 1.96 * sd(estimate)
    )
  d1$x <- 1:nrow(d1)
  
  p1data <- mySSF %>% 
    mutate(coef = map(model, ~ broom::tidy(.x$model))) %>%
    select(ID, Sex, uphillRate, coef, mvmtClass) %>% 
    unnest(cols = c(coef)) %>%
    mutate(ID = factor(ID)) %>%
    mutate(termFac = as.factor(term))
  
  otlrs = p1data %>%
    filter(estimate < -10 |
             estimate > 10)
  
  if(nrow(otlrs)>0){
    message("Outlier estimates:")
    print(otlrs)
  }
  
  p1 = p1data  %>%
    arrange(uphillRate) %>%
    ggplot(aes(x = termFac, y = exp(estimate))) +
    geom_point(col = "transparent") +
    geom_hline(yintercept = 1, 
               lty = 2) +
    geom_point(aes(col = uphillRate,
                   group = ID), 
               position = position_dodge(width = 0.5)) +
    scale_color_viridis_c("Rate of migration") +
    labs(x = "Habitat", y = "Relative Selection Strength") +
    scale_x_discrete() +
    CJsBasics::BasicTheme +
    facet_wrap(~Sex) +
    theme(axis.text.x = element_text(angle = 90,vjust=0.5,hjust=1))
  if(drawplot==TRUE){plot(p1)}
  return(p1data)
}
## Plotting migration rate
plot_rateInteractions = function(mySSF){
  coeffs = plot_popCoeffs(mySSF,drawplot=FALSE)
  p2 = coeffs %>%
    ggplot(aes(x = uphillRate, y = estimate)) +
    geom_hline(aes(yintercept = 0),
               lty = 2,
               col = "grey60") +
    geom_errorbar(aes(x = uphillRate,
                      ymin = estimate-std.error,
                      ymax = estimate+std.error,
                      group = ID)) +
    geom_point(aes(col = Sex)) +
    stat_smooth(col = "grey40",
                method = "lm") +
    facet_wrap(~term) +
    labs(x = "Rate of migration",
         y = "Relative selection strength") +
    CJsBasics::BasicTheme
  plot(p2)
}

#### Import data ----
## Extracts file paths
allExtracts_FPs = list.files(path = "data/extract",
                             pattern = "cov_extract_*",
                             full.names = TRUE)
## Import extracts
allExtracts_list = lapply(allExtracts_FPs,
                          function(x){
                            y = read_csv(x,
                                         col_types = cols(Sex = "c"))
                            return(y)
                          })
for(i in 1:length(allExtracts_list)){
  print(unique(allExtracts_list[[i]]$Sex))
}
## Migration patterns 
migPhenol = read_csv("data/mig_PatternsCompiled.csv") %>%
  rename(Individual = ID) %>%
  mutate(ID = paste0(Individual,"_",springYr))
indParams = migPhenol %>%
  dplyr::select(ID, Sex, mvmtClass, upDOY, uphillRate)

#### Compile the extract data ----
myDat = do.call("rbind",
                allExtracts_list)

#### Check for 0-variance variables ----
## This pulls out just endpoint columns for each step and summarizes its variance
myVars = myDat %>%
  nest(-stepID) %>%
  mutate(var = lapply(data,
                      function(x){
                        xnames = names(x)
                        colsStart = grep("_start", xnames)
                        colsEnd = grep("_end", xnames)
                        endVals = x[,colsEnd] %>%
                          select(-tod_end_)
                        endVals$ID = x$ID
                        endVals$Year = x$Year
                        
                        endVar = endVals %>%
                          pivot_longer(-c(ID,Year)) %>%
                          group_by(ID,Year,name) %>%
                          summarize(var = var(value, na.rm = T)) %>%
                          add_column(point = "end") %>%
                          ungroup()
                        
                        return(endVar)
                      }))

## Check for individual/variable combinations where the 
# endpoint variance is 0 in greater than 90% of steps
myVars %>%
  dplyr::select(stepID, var) %>%
  unnest(var) %>%
  group_by(name, ID) %>%
  dplyr::summarize(lowVarRate = sum(var==0,
                                    na.rm = T)/n()) %>%
  dplyr::filter(lowVarRate > 0.80) %>%
  pull(name) %>%
  unique()

## Now, remove cases with snow phenology flags, and scale/center predictor variables.
myDat_wrangled_mostly = myDat %>%
  mutate(
    ## Snow variables
    ## Remove the flagged snow phenology points
    meltMidpoint_start = case_when(meltMidpoint_start > 996 ~ as.numeric(NA),
                                   TRUE ~ meltMidpoint_start),
    meltMidpoint_end = case_when(meltMidpoint_end > 996 ~ as.numeric(NA),
                                 TRUE ~ meltMidpoint_end),
    fallMidpoint_start = case_when(fallMidpoint_start > 996 ~ as.numeric(NA),
                                   TRUE ~ fallMidpoint_start),
    fallMidpoint_end = case_when(fallMidpoint_end > 996 ~ as.numeric(NA),
                                 TRUE ~ fallMidpoint_end),
    # Fractional snow cover rescaled to [0,1] rather than [0,100]
    # this makes more sense than centering/SD scaling because 0 = 0 snow
    FSC_end_scaled = FSC_end/100,
    # Snow conditions: how far traversed this step?
    deltaSnowDist = snowDist_end - snowDist_start,
    # Snow conditions at start and end of step
    snowDist_start_scaled = log(snowDist_start + 1),
    snowDist_end_scaled = log(snowDist_end + 1),
    # Date of snowfall at start and end of step
    snowFall_start = lubridate::ymd(paste0(Year-1, "-", startMonthDay)) + 
      lubridate::days(round(fallMidpoint_start) - 1),
    snowFall_end = lubridate::ymd(paste0(Year-1, "-", startMonthDay)) + 
      lubridate::days(round(fallMidpoint_end) - 1),
    # Date of snowmelt at start and end of step
    snowMelt_start = lubridate::ymd(paste0(Year-1, "-", startMonthDay)) + 
      lubridate::days(round(meltMidpoint_start) - 1),
    snowMelt_end = lubridate::ymd(paste0(Year-1, "-", startMonthDay)) + 
      lubridate::days(round(meltMidpoint_end) - 1),
    # Difference between GPS loc time and snowmelt timing
    deltaMelt_start = as.numeric(difftime(as.Date(t1_),
                                          as.Date(snowMelt_start),
                                          units = "days")),
    deltaMelt_end = as.numeric(difftime(as.Date(t1_),
                                        as.Date(snowMelt_end),
                                        units = "days")),
    
    ## Forage conditions
    # Green-up timing
    greenUp_start = lubridate::ymd(paste0(Year-1, "-", startMonthDay)) + 
      lubridate::days(round(greenupMidpoint_start) - 1),
    greenUp_end = lubridate::ymd(paste0(Year-1, "-", startMonthDay)) + 
      lubridate::days(round(greenupMidpoint_end) - 1),
    # Difference between GPS loc time and green-up timing
    deltaGreen_start = as.numeric(difftime(as.Date(t1_),
                                           as.Date(greenUp_start),
                                           units = "days")),
    deltaGreen_end = as.numeric(difftime(as.Date(t1_),
                                         as.Date(greenUp_end),
                                         units = "days")),
    
    # Scaled NDVI
    NDVI_start_scaled = ndviCurve(DOY = as.numeric(strftime(t1_,
                                                            format = "%j")),
                                  xmidS = greenupMidpoint_start,
                                  xmidA = senescMidpoint_start,
                                  scalS = greenupRate_start,
                                  scalA = senescRate_start),
    NDVI_end_scaled = ndviCurve(DOY = as.numeric(strftime(t1_,
                                                          format = "%j")),
                                xmidS = greenupMidpoint_end,
                                xmidA = senescMidpoint_end,
                                scalS = greenupRate_end,
                                scalA = senescRate_end),
    
    ## Weather conditions
    highTemp_start_scaled = scale(highTemp_start),
    highTemp_end_scaled = scale(highTemp_end),
    highTemp_start_over27 = as.factor(ifelse(highTemp_start > 27,
                                             "Yes","No")),
    
    Precip_start_scaled = log(Precip_start + 1),
    Precip_end_scaled = log(Precip_end + 1),
    raining_start = as.factor(ifelse(Precip_start > 0,
                                     "Yes","No")),
    is_raining_start_numeric = ifelse(raining_start == "Yes",
                                      1,
                                      0),
    ## Terrain conditions
    # Movement across topography
    deltaDEM = DEM_end - DEM_start,
    deltaDEM_scaled = scale(deltaDEM),
    # Scaled elevation at start and endpoints
    DEM_start_scaled = scale(DEM_start),
    DEM_end_scaled = scale(DEM_end),
    # Scaled slope at start and endpoints
    SLP_start_scaled = scale(SLP_start),
    SLP_end_scaled = scale(SLP_end),
    # Cosine of aspect at start and endpoints
    cos_ASP_start = cos((ASP_start*pi)/180),
    cos_ASP_end = cos((ASP_end*pi)/180),
    # Distance from escape terrain at start and endpoints
    ESC_Dist_start_scaled = log(ESC_Dist_start + 1),
    ESC_Dist_end_scaled = log(ESC_Dist_end + 1),
    
    ## Land cover conditions
    NLCD_higherOrder_start = substr(NLCD_start, start = 1, stop = 1),
    NLCD_higherOrder_end = substr(NLCD_end, start = 1, stop = 1),
    NLCD_start = as.factor(NLCD_start),
    NLCD_end = as.factor(NLCD_end),
    
    ## Movement parameters
    cos_ta_ = cos(ta_),
    log_sl_ = log(sl_),
    case_ = as.numeric(case_)) %>%
  group_by(ID, Year) %>%
  mutate(is_hotDay_start = highTemp_start > quantile(highTemp_start, 0.75),
         is_hotDay_start_numeric = as.numeric(is_hotDay_start),
         is_hotDay_start_factor = as.character(is_hotDay_start_numeric)) %>%
  ungroup()

# Scales doesn't like vectorization I guess so use rowwise for this
# Computers hate it!
myDat_wrangled = myDat_wrangled_mostly %>%
  rowwise() %>%
  mutate(
    # Absolute NDVI
    NDVI_start_absol = ifelse(!is.na(NDVI_start_scaled),
                              scales::rescale(NDVI_start_scaled,
                                              from = c(0,1),
                                              to = c(wNDVI_start,
                                                     mNDVI_start)),
                              as.numeric(NA))/10000,
    NDVI_end_absol = ifelse(!is.na(NDVI_end_scaled), 
                            scales::rescale(NDVI_end_scaled,
                                            from = c(0,1),
                                            to = c(wNDVI_end,
                                                   mNDVI_end)),
                            as.numeric(NA))/10000
  ) %>%
  ungroup()

#### Save dataset ----
saveRDS(object = myDat_wrangled,
        file = "data/ssfData.RDS")


#### Epilogue ----
stopTime = Sys.time()
message("elapsed time: ")
stopTime-startTime
BRRR::skrrrahh(12)
Sys.sleep(5)
## Software versions
sessionInfo()

