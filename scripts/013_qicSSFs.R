#### Load packages ----
library(tidyverse)
library(amt)
library(survival)

#### Import data ----
## Extracted step selection features
## Use only complete cases because model comparison relies on consistent observations
myDat = readRDS("data/ssfData.RDS") %>%
  filter(complete.cases(.)) %>%
  mutate(Year = as.factor(Year))


#### Fit step selection functions ----
## Comparison by QIC (a quasi-likelihood weighted generalized measure for AIC [20]) among models fit with all predictors and nested models fit with only terrain, weather, snow, or forage parameters revealed the most support for the overall model (∆QIC > 500 for all nested models). 
## Create an overall model and then terrain, weather, snow, and forage models
ssf1 = clogit(case_ ~ 
                DEM_end_scaled:is_hotDay_start_numeric +
                DEM_end_scaled:is_raining_start_numeric +
                DEM_end_scaled +
                SLP_end_scaled +
                cos_ASP_end + 
                ESC_Dist_end_scaled +
                highTemp_end_scaled +
                Precip_end_scaled +
                FSC_end_scaled + 
                snowDist_end_scaled +
                deltaMelt_end +
                NDVI_end_absol +
                NDVI_end_scaled +
                deltaGreen_end +
                strata(stepID),
              data = myDat,
              model = TRUE,
              robust = TRUE,method = "efron")
ssf1terrain = clogit(case_ ~ 
                       DEM_end_scaled +
                       SLP_end_scaled +
                       cos_ASP_end + 
                       ESC_Dist_end_scaled +
                       strata(stepID),
                     data = myDat,
                     model = TRUE,
                     robust = TRUE,method = "efron")
ssf1weather = clogit(case_ ~ 
                       highTemp_end_scaled +
                       Precip_end_scaled +
                       strata(stepID),
                     data = myDat,
                     model = TRUE,
                     robust = TRUE,method = "efron")
ssf1snow = clogit(case_ ~ 
                    FSC_end_scaled + 
                    snowDist_end_scaled +
                    deltaMelt_end +
                    strata(stepID),
                  data = myDat,
                  model = TRUE,
                  robust = TRUE,method = "efron")
ssf1forage = clogit(case_ ~ 
                      NDVI_end_absol +
                      NDVI_end_scaled +
                      deltaGreen_end +
                      strata(stepID),
                    data = myDat,
                    model = TRUE,
                    robust = TRUE,method = "efron")


#### Compare models ----
lmtest::lrtest(ssf1,
               ssf1terrain)

# Model comparison: From Justine's site
QIC.coxph <- function(mod, details = FALSE) {
  trace <- sum(diag(solve(mod$naive.var) %*% mod$var))
  quasi <- mod$loglik[2]
  return(-2*quasi + 2*trace)
}
# Which model is best supported?
myQICs = c(QIC.coxph(ssf1),
           QIC.coxph(ssf1terrain),
           QIC.coxph(ssf1weather),
           QIC.coxph(ssf1snow),
           QIC.coxph(ssf1forage))
myQICs - myQICs[1]

myAICs = MuMIn::AICc(ssf1, ssf1terrain, ssf1weather, ssf1snow, ssf1forage)$AICc
myAICs - myAICs[1]



modelList = list(ssf1, ssf1terrain, ssf1weather, ssf1snow, ssf1forage)
modelNames = c("complete", "terrain", 
               "weather", "snow", "forage")

generateMultiSSFdata = function(mySSF){
  # Visualize variation in selection propensity at the population level
  p1data <- as.data.frame(summary(mySSF)$coefficients) %>%
    mutate("term" = rownames(.),
           log_coef = exp(coef),
           log_min = exp(coef - `se(coef)`),
           log_max = exp(coef + `se(coef)`)) %>%
    mutate(termRenamed = case_when(term == "DEM_end_scaled:is_hotDay_start_numeric" ~
                                     "Elevation*Hot day (start)",
                                   term == "DEM_end_scaled:hot_start" ~
                                     "Elevation * Hot (start)",
                                   term == "DEM_end_scaled:highTemp_start" ~
                                     "Elevation * Temp (start)",
                                   term == "DEM_end_scaled:highTemp_start_scaled" ~
                                     "Elevation * Temp (start)",
                                   term == "DEM_end_scaled:is_raining_start_numeric" ~
                                     "Elevation*Raining (start)",
                                   term == "DEM_end_scaled:Precip_start" ~
                                     "Elevation * Precip (start)",
                                   term == "DEM_end_scaled:Precip_start_scaled" ~
                                     "Elevation * Precip (start)",
                                   term == "DEM_end_scaled:FSC_start" ~
                                     "Elevation * FSC (start)",
                                   term == "DEM_end_scaled:NDVI_end_scaled" ~
                                     "Elevation * NDVI",
                                   term == "DEM_end_scaled" ~
                                     "Elevation",
                                   term == "highTemp_end_scaled" ~
                                     "Temperature",
                                   term == "highTemp_start_scaled" ~
                                     "Temperature (start)",
                                   term == "Precip_end_scaled" ~
                                     "Precipitation",
                                   term == "Precip_start_scaled" ~
                                     "Precipitation (start)",
                                   term == "snowDist_end_scaled" ~
                                     "Distance to snow",
                                   term == "FSC_end_scaled" ~
                                     "Fractional snow cover",
                                   term == "FSC_end" ~
                                     "Fractional snow cover",
                                   term == "deltaMelt_end" ~
                                     "Days since peak snowmelt",
                                   term == "deltaMelt_end_months" ~
                                     "Months since peak FSC depletion",
                                   term == "deltaMelt_end_abs" ~
                                     "Time from peak snowmelt",
                                   term == "deltaMelt_end_months_abs" ~
                                     "Time from peak snowmelt",
                                   term == "I(deltaMelt_end^2)" ~
                                     "(Days since peak snowmelt)^2",
                                   term == "NDVI_end_scaled" ~
                                     "Relative NDVI",
                                   term == "NDVI_end_absol" ~
                                     "Absolute NDVI",
                                   term == "deltaGreen_end" ~
                                     "Days since peak green-up",
                                   term == "deltaGreen_end_months" ~
                                     "Months since peak green-up",
                                   term == "deltaGreen_end_abs" ~
                                     "Time from peak green-up",
                                   term == "deltaGreen_end_months_abs" ~
                                     "Time from peak green-up",
                                   term == "I(deltaGreen_end^2)" ~
                                     "(Days since peak green-up)^2",
                                   term == "DEM_end_scaled" ~
                                     "Elevation",
                                   term == "SLP_end_scaled" ~
                                     "Terrain slope",
                                   term == "cos_ASP_end" ~
                                     "Terrain northness",
                                   term == "ESC_Dist_end_scaled" ~
                                     "Distance to escape terrain",
                                   term == "sl_" ~
                                     "Step length",
                                   term == "cos_ta_" ~
                                     "Turning angle (cos-transformed)",
                                   TRUE ~ as.character(NA)
    )) %>%
    mutate(term = factor(term,
                         levels = c(
                           ## Movement drivers: temp, precip, snow, plants
                           "DEM_end_scaled:is_hotDay_start_numeric",
                           "DEM_end_scaled:highTemp_start",
                           "DEM_end_scaled:highTemp_start_scaled",
                           "DEM_end_scaled:raining_start",
                           "DEM_end_scaled:Precip_start",
                           "DEM_end_scaled:Precip_start_scaled",
                           "DEM_end_scaled:is_raining_start_numeric",
                           "DEM_end_scaled:FSC_start",
                           "DEM_end_scaled:deltaMelt_start",
                           "DEM_end_scaled:NDVI_start_scaled",
                           "DEM_end_scaled:deltaGreen_start",
                           "DEM_end_scaled:NDVI_end_scaled",
                           # Selection factors: temp, snow, plants, terrain
                           # - Temp
                           "highTemp_end_scaled",
                           "highTemp_start_scaled",
                           # - Precip
                           "Precip_end_scaled",
                           "Precip_start_scaled",
                           # - Snow
                           "snowDist_end_scaled",
                           "FSC_end",
                           "FSC_end_scaled",
                           "deltaMelt_end",
                           "deltaMelt_end_abs",
                           "deltaMelt_end_months",
                           "deltaMelt_end_months_abs",
                           # - Plants
                           "NDVI_end_scaled",
                           "NDVI_end_absol",
                           "deltaGreen_end",
                           "deltaGreen_end_abs",
                           "deltaGreen_end_months",
                           "deltaGreen_end_months_abs",
                           # - Terrain
                           "DEM_end_scaled", 
                           "SLP_end_scaled",
                           "cos_ASP_end", 
                           "ESC_Dist_end_scaled",
                           "sl_", "log_sl_", "cos_ta_")))  %>%
    mutate(termFac = factor(termRenamed,
                            levels = c("Elevation*Hot day (start)",
                                       "Elevation * Temp (start)",
                                       "Elevation*Raining (start)",
                                       "Elevation * Precip (start)",
                                       "Elevation * FSC (start)",
                                       "Elevation * NDVI",
                                       "Elevation",
                                       "Terrain slope",
                                       "Terrain northness",
                                       "Distance to escape terrain",
                                       "Temperature",
                                       "Temperature (start)",
                                       "Precipitation",
                                       "Precipitation (start)",
                                       "Fractional snow cover",
                                       "Distance to snow",
                                       "Time from peak snowmelt",
                                       "Days since peak snowmelt",
                                       "Months since peak FSC depletion",
                                       "(Days since peak snowmelt)^2",
                                       "Absolute NDVI",
                                       "Relative NDVI",
                                       "Time from peak green-up",
                                       "Days since peak green-up",
                                       "Months since peak green-up",
                                       "(Days since peak green-up)^2",
                                       "Step length",
                                       "Turning angle (cos-transformed)")))
  p1data
}

multissfdatalist = list()
for(i in 1:length(modelList)){
  multissfdatalist[[i]] = generateMultiSSFdata(modelList[[i]]) %>%
    mutate(model = modelNames[i])
}

multiSSF_data = do.call("rbind", multissfdatalist)

multiSSF_data   %>%
  mutate(model = factor(model, levels = modelNames),
         sigpos = ifelse((log_min > 1),
                         "YES",
                         "NO"),
         signeg = ifelse((log_max < 1),
                         "YES",
                         "NO"),
         sigerr = ifelse(sigpos == "YES", 1, ifelse(signeg == "YES", 1, 0.5))) %>% 
  arrange(desc(model)) %>%
  ggplot() +
  geom_hline(aes(yintercept=1),
             lty=2,col="grey40") +
  geom_errorbar(aes(x = termFac,
                    ymin = log_min,
                    ymax = log_max,
                    group = model),
                col = "black",
                position = position_dodge(width = 0.5),
                size = 2,
                width = 0.75) +
  geom_errorbar(aes(x = termFac,
                    ymin = log_min,
                    ymax = log_max,
                    group = model),
                col = "white",
                position = position_dodge(width = 0.5),
                size = 1,
                width = 0.75) +
  geom_errorbar(aes(x = termFac,
                    ymin = log_min,
                    ymax = log_max,
                    col = model),
                alpha = 0.75,
                position = position_dodge(width = 0.5),
                size = 1,
                width = 0.75) +
  geom_point(aes(x = termFac,
                 y = log_coef,
                 col = model),
             position = position_dodge(width = 0.5),
             size = 0.75) +
  scale_color_manual("Model",
                     values = c("black",
                                "brown","yellow3",
                                "lightblue","forestgreen")) +
  xlab("") +
  ylab(expression(paste("Model ", beta))) +
  CJsBasics::BasicTheme +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1, 
                                   vjust = 0.5),
        legend.key.height =  unit(0.01,"in"), 
        legend.position = c(0.15,0.15))
