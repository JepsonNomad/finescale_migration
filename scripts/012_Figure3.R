#### Load packages ----
library(tidyverse)
library(cowplot)
library(amt)
library(kableExtra)
library(hab)
library(survival)


#### Define functions ----
source("scripts/000_ssfPlottingFx.R")


#### Import data ----
mig = read_csv("data/mig_PatternsCompiled.csv")

## Extracted step selection features
myDat_raw = readRDS("data/ssfData.RDS")

## Note that hab::kfold fails if there are any NA's
## So we need to remove groups that have incomplete cases.
## There has got to be a more elegant way of coding this.
myDat = myDat_raw %>% 
  dplyr::select(ID,
                Sex,
                HU,
                Year,
                MvmtClass,
                burst_,
                case_,
                tod_start_,
                tod_end_,
                stepID,
                is_hotDay_start_numeric ,
                is_raining_start_numeric ,
                DEM_end_scaled ,
                SLP_end_scaled ,
                cos_ASP_end , 
                ESC_Dist_end_scaled ,
                highTemp_end_scaled ,
                Precip_end_scaled ,
                FSC_end_scaled , 
                snowDist_end_scaled ,
                deltaMelt_end,
                NDVI_end_absol ,
                NDVI_end_scaled ,
                deltaGreen_end) %>%
  dplyr::group_by(stepID) %>% 
  filter(!any(is.na(c(ID,
                      Sex,
                      HU,
                      Year,
                      MvmtClass,
                      burst_,
                      case_,
                      tod_start_,
                      tod_end_,
                      stepID,
                      is_hotDay_start_numeric,
                      is_raining_start_numeric,
                      DEM_end_scaled,
                      SLP_end_scaled,
                      cos_ASP_end,
                      ESC_Dist_end_scaled,
                      highTemp_end_scaled,
                      Precip_end_scaled,
                      FSC_end_scaled,
                      snowDist_end_scaled,
                      deltaMelt_end,
                      NDVI_end_absol,
                      NDVI_end_scaled,
                      deltaGreen_end)))) %>%
  ungroup()
myDat

#### Fit step selection function ----
ssf1 = myDat %>%
  mutate(Year = as.factor(Year)) %>%
  fit_clogit(case_ ~ 
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
             model = TRUE)
summary(ssf1)
class(ssf1)
ssf1_clogit = as(ssf1[[1]], "clogit")
ssf1_kfold = hab::kfold(ssf1_clogit, 
                        reproducible = T,
                        jitter = T,
                        k = 5,
                        nrepet = 100)
ssf1_kfold
ssf1_kfold %>%
  filter(type == "obs") %>%
  pull(cor) %>%
  mean()
# [1] 0.9340486
## with sd instead of mean:
# [1] 0.01909392

ssf1_kfold %>%
  filter(type != "obs") %>%
  pull(cor) %>%
  mean()
# [1] -0.01541244
## with sd instead of mean:
# [1] 0.1873906

ssf1_kfold_plot = ssf1_kfold %>% 
  mutate(type = ifelse(type == "obs",
                       "Observed",
                       "Random")) %>%
  ggplot() +
  geom_hline(aes(yintercept=0),col="grey40",lty=2) +
  geom_boxplot(aes(x = type, y = cor),
               fill = "white") +
  xlab("Group") +
  ylab("Correlation") +
  scale_y_continuous(limits = c(-1,1),
                     breaks=c(-1,0,1),
                     labels=c("-1","0","1")) +
  CJsBasics::BasicTheme
ggsave("plots/ssf_kfold_correlations.jpg", ssf1_kfold_plot,
       width = 2.5, height = 4, units = "in", dpi = 300)

allBighorn_fullModel = plotSSF(ssf1) + 
  ylab(expression(paste("Model ", beta))) 
allBighorn_fullModel


#### Plotting RSS ----
## Generate plot
noIntPltSSF = plotSSF(ssf1) + 
  ylab(expression(paste("Model ", beta))) +
  theme(plot.background = element_rect(fill="transparent"),
        panel.background = element_rect(fill="transparent"))
noIntPltRSS1 = plotRSS_snow(ssf1, md = myDat) + 
  theme(plot.background = element_rect(fill="transparent"),
        panel.background = element_rect(fill="transparent"))
noIntPltRSS2 = plotRSS_green(ssf1, md = myDat) + 
  theme(plot.background = element_rect(fill="transparent"),
        panel.background = element_rect(fill="transparent"))

noIntPltSSF
noIntPltRSS1
noIntPltRSS2

compiled_RSS = plot_grid(noIntPltSSF,
                         plot_grid(noIntPltRSS1 + 
                                     theme(
                                       plot.margin = margin(t = 0.1,
                                                            r = 0.05,
                                                            l = 0.25,
                                                            unit = "in")),
                                   noIntPltRSS2 + 
                                     theme(
                                       plot.margin = margin(r = 0.05,
                                                            b = 0.1,
                                                            l = 0.25,
                                                            unit = "in")),
                                   nrow = 2,
                                   align = "hv",
                                   labels = c("B","C")),
                         labels = c("A",""),
                         nrow = 1,
                         rel_widths = c(0.85,1))
compiled_RSS
ggsave("plots/compiledSelection.jpg",
       compiled_RSS,
       width = 6, height = 4.5, units = "in", dpi = 300)
ggsave("plots/compiledSelection.svg",
       compiled_RSS,
       width = 6, height = 4.5, units = "in", dpi = 300)


#### Visualize the effect of a storm ----
## Supplementary materials S4
rainMetadata = myDat %>%
  filter(is_raining_start_numeric > 0) %>%
  filter(Precip_end_scaled > 2) %>%
  pull(stepID) %>%
  unique()
rainInds = sapply(str_split(rainMetadata, pattern = "_"),
                  "[[", 1)
rainYrs = sapply(str_split(rainMetadata, pattern = "_"),
                 "[[", 2)

rainData = myDat_raw %>%
  mutate(indYrFinder = paste0(ID, "_", Year)) %>%
  filter(indYrFinder %in% paste0(rainInds, "_", rainYrs))

rainDates = rainData %>%
  filter(case_ == 1,
         Precip_start > 0) %>%
  group_by(HU, Year) %>%
  summarize(dates = paste0(t1_, collapse = ", "))

huSelections = rainData %>%
  mutate(rainyHUyrs = paste0(rainData$HU, " ",
                             rainData$Year)) %>%
  pull(rainyHUyrs) %>%
  unique()
huPlotter = function(myHU = NULL,
                     myYear = NULL){
  xintsHMS = rainDates %>% 
    filter(HU == myHU,
           Year == myYear) %>%
    pull(dates) %>%
    str_split(pattern = ", ") %>%
    unlist() %>%
    lubridate::ymd_hms()
  xints = as.Date(xintsHMS) %>%
    unique() %>%
    paste0(" 00:00:01 UTC") %>%
    lubridate::ymd_hms()
  # xints
  rainElevPlot = rainData %>%
    rename(Individual = ID) %>%
    filter(HU == myHU,
           Year == myYear,
           case_ == 1) %>%
    ggplot(aes(x = lubridate::ymd_hms(t1_),
               y = DEM_end,
               group = Individual)) +
    geom_vline(xintercept = xints,
               lwd = 1, col = "grey80") +
    geom_line(aes(linetype = Individual)) +
    ggtitle(paste0(myHU, " ", myYear)) +
    xlab("Date") +
    ylab("Elevation (m)") +
    CJsBasics::BasicTheme +
    theme(legend.position = c(0.9,0.2),
          legend.key.height = unit(0.15,"in"),
          legend.background = element_blank())
  return(rainElevPlot)
}
for(i in huSelections){
  huI = str_split(i,pattern = " ")[[1]][1]
  yrI = str_split(i,pattern = " ")[[1]][2]
  plot(huPlotter(myHU = huI, myYear = yrI))
}
exampleRainDrivenMvmt = huPlotter(myHU = "Cn",
                                  myYear = 2016)
exampleRainDrivenMvmt
ggsave(filename = "plots/figure_S4.jpg",
       plot = exampleRainDrivenMvmt,
       width = 4.5, height = 3, units = "in",dpi = 300)