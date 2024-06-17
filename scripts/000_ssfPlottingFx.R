
plotSSF = function(mySSF){
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
                                     "Months since peak snowmelt",
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
                           "DEM_end_scaled:hot_start",
                           "DEM_end_scaled:highTemp_start",
                           "DEM_end_scaled:highTemp_start_scaled",
                           "DEM_end_scaled:is_raining_start_numeric",
                           "DEM_end_scaled:Precip_start",
                           "DEM_end_scaled:Precip_start_scaled",
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
                           "NDVI_end_absol",
                           "NDVI_end_scaled",
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
                                       "Elevation*Not hot day (start)",
                                       "Elevation * Hot (start)",
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
                                       "Relative NDVI",
                                       "Absolute NDVI",
                                       "Time from peak green-up",
                                       "Days since peak green-up",
                                       "Months since peak green-up",
                                       "(Days since peak green-up)^2",
                                       "Step length",
                                       "Turning angle (cos-transformed)")))
  
  p1 = p1data  %>%
    mutate(sigpos = ifelse((log_min > 1),
                           "YES",
                           "NO"),
           signeg = ifelse((log_max < 1),
                           "YES",
                           "NO"),
           sigerr = ifelse(sigpos == "YES", 1, ifelse(signeg == "YES", 1, 0.5))) %>% 
    ggplot() +
    geom_hline(aes(yintercept=1),lty=2,col="grey40") +
    geom_errorbar(aes(x = termFac,
                      ymin = log_min,
                      ymax = log_max,
                      alpha = sigerr),
                  width = 0.5) +
    geom_point(aes(x = termFac,
                   y = log_coef,
                   alpha = sigerr),
               size = 0.75) +
    scale_alpha_identity() +
    # geom_errorbar(aes(x = termFac,
    #                   ymin = log_min,
    #                   ymax = log_max,
    #                   col = sigpos),
    #               lwd = 2,
    #               width = 0,
    #               alpha = 0.5,
    #               show.legend = F) +
    # scale_color_manual(values = c("transparent", "royalblue3")) +
    # ggnewscale::new_scale_color() +
    # geom_errorbar(aes(x = termFac,
  #                   ymin = log_min,
  #                   ymax = log_max,
  #                   col = signeg),
  #               lwd = 2,
  #               width = 0,
  #               alpha = 0.5,
  #               show.legend = F) +
  # scale_color_manual(values = c("transparent", "tomato2")) +
  xlab("") +
    ylab("Relative selection strength") +
    CJsBasics::BasicTheme +
    theme(axis.text.x = element_text(angle = 90,
                                     hjust = 1, 
                                     vjust = 0.5))
  
  
  p1
}
plotRSS_snow = function(mySSF, md){
  x1 <- data.frame(is_hotDay_start_numeric = 0,
                   is_raining_start_numeric = 0,
                   highTemp_end_scaled =0,
                   snowDist_end_scaled =0,
                   FSC_end =0, 
                   FSC_end_scaled =0, 
                   Precip_end_scaled = 0,
                   deltaMelt_end = seq(from = 0,
                                       to = 200,
                                       length.out = 500),
                   deltaMelt_end_abs = seq(from = 0,
                                       to = 200,
                                       length.out = 500),
                   deltaMelt_end_months_abs = seq(from = 0,
                                                  to = 8,
                                                  length.out = 500),
                   NDVI_end_absol = 0,
                   NDVI_end_scaled =0,
                   deltaGreen_end =0,
                   deltaGreen_end_abs =0,
                   `abs(deltaGreen_end)` =0,
                   deltaGreen_end_months_abs = 0,
                   DEM_end_scaled =0,
                   SLP_end_scaled =0,
                   cos_ASP_end =0, 
                   ESC_Dist_end_scaled =0,
                   sl_ = 0,
                   cos_ta_ = 0,
                   stepID = "S18_2004_9833")
  x1$`I(deltaMelt_end^2)` = I(x1$deltaMelt_end)^2
  
  
  x2 <- data.frame(is_hotDay_start_numeric = 0,
                   is_raining_start_numeric = 0,
                   highTemp_end_scaled = 0,
                   snowDist_end_scaled = 0,
                   FSC_end = 0, 
                   FSC_end_scaled =0, 
                   Precip_end_scaled = 0,
                   deltaMelt_end = 0,
                   deltaMelt_end_abs = 0,
                   deltaMelt_end_months_abs = 0,
                   NDVI_end_absol = 0,
                   NDVI_end_scaled = 0,
                   deltaGreen_end = 0,
                   deltaGreen_end_abs = 0,
                   `abs(deltaGreen_end)` = 0,
                   deltaGreen_end_months_abs = 0,
                   DEM_end_scaled = 0,
                   SLP_end_scaled = 0,
                   cos_ASP_end = 0, 
                   ESC_Dist_end_scaled =0,
                   sl_ = 0,
                   cos_ta_ = 0,
                   stepID = "S18_2004_9833")
  x2$`abs(deltaMelt_end)` = mean(I(abs(x2$deltaMelt_end)), na.rm = T)
  
  logRSS <- log_rss(mySSF, x1, x2, ci = "se", ci_level = 0.95)
  
  # We have a plot method for 'log_rss' objects to make a very basic figure.
  # plot(logRSS)
  
  # But if we want more control, we can use ggplot with the 'df' data.frame.
  ggplot(logRSS$df, aes(x = deltaMelt_end_x1,
                        y = exp(log_rss), 
                        ymin = exp(lwr), 
                        ymax = exp(upr))) +
    geom_ribbon(color = "transparent", fill = "grey80") +
    geom_line(color = "#00bcfc") +
    geom_hline(yintercept = 1, color = "grey40", linetype = "dashed") +
    xlab("Days since peak snowmelt") +
    ylab("Selection strength") +
    CJsBasics::BasicTheme
}
plotRSS_green = function(mySSF, md){
  x1 <- data.frame(is_hotDay_start_numeric = 0,
                   is_raining_start_numeric = 0,
                   highTemp_end_scaled = 0,
                   snowDist_end_scaled = 0,
                   FSC_end =0, 
                   FSC_end_scaled =0, 
                   Precip_end_scaled = 0,
                   deltaMelt_end = 0,
                   deltaMelt_end_abs = 0,
                   `abs(deltaMelt_end)` = 0,
                   deltaMelt_end_months_abs = 0,
                   NDVI_end_absol = 0,
                   NDVI_end_scaled = 0,
                   deltaGreen_end = seq(from = 0,
                                        to = 200,
                                        length.out = 500),
                   deltaGreen_end_abs = seq(from = 0,
                                        to = 200,
                                        length.out = 500),
                   deltaGreen_end_months_abs = seq(from = 0,
                                                   to = 8,
                                                   length.out = 500),
                   DEM_end_scaled =0,
                   SLP_end_scaled =0,
                   cos_ASP_end =0, 
                   ESC_Dist_end_scaled =0,
                   sl_ = 0,
                   cos_ta_ = 0,
                   stepID = "S18_2004_9833")
  x1$`I(deltaGreen_end^2)` = I(x1$deltaGreen_end)^2
  
  
  x2 <- data.frame(is_hotDay_start_numeric = 0,
                   is_raining_start_numeric = 0,
                   highTemp_end_scaled = 0,
                   snowDist_end_scaled = 0,
                   FSC_end = 0, 
                   FSC_end_scaled =0, 
                   Precip_end_scaled = 0,
                   deltaMelt_end = 0,
                   deltaMelt_end_abs = 0,
                   `abs(deltaMelt)` = 0,
                   deltaMelt_end_months_abs = 0,
                   NDVI_end_absol = 0,
                   NDVI_end_scaled = 0,
                   deltaGreen_end = 0,
                   deltaGreen_end_abs = 0,
                   deltaGreen_end_months_abs = 0,
                   DEM_end_scaled = 0,
                   SLP_end_scaled = 0,
                   cos_ASP_end = 0, 
                   ESC_Dist_end_scaled = 0,
                   sl_ = 0,
                   cos_ta_ = 0,
                   stepID = "S18_2004_9833")
  x2$`abs(deltaGreen_end)` = mean(I(abs(x2$deltaGreen_end)), na.rm = T)
  
  logRSS <- log_rss(mySSF, x1, x2, ci = "se", ci_level = 0.95)
  
  # We have a plot method for 'log_rss' objects to make a very basic figure.
  # plot(logRSS)
  
  # But if we want more control, we can use ggplot with the 'df' data.frame.
  ggplot(logRSS$df, aes(x = deltaGreen_end_x1,
                        y = exp(log_rss), 
                        ymin = exp(lwr), 
                        ymax = exp(upr))) +
    geom_ribbon(color = "transparent", fill = "grey80") +
    geom_line(color = "#00b100") +
    geom_hline(yintercept = 1, color = "grey40", linetype = "dashed") +
    xlab("Days since peak green-up") +
    ylab("Selection strength") +
    CJsBasics::BasicTheme
}

