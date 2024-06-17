#### Load packages ----
library(tictoc)
library(tidyverse)
library(cowplot)
library(amt)
library(kableExtra)

#### Define functions ----
source("scripts/000_ssfPlottingFx.R")

#### Import data ----
## Extracted step selection features
myDat = readRDS("data/ssfData.RDS")

migPat = read_csv("data/mig_PatternsCompiled.csv") %>%
  select(ID, Sex, springYr, mvmtClass, upDOY, uphillRate, highElev, delta)

#### What do different types of movements accomplish? ----
#### Migrant vs resident movements ----
## Migrants
migrIndYrs = myDat %>%
  left_join(migPat %>% mutate(Year = springYr),
            by = c("ID", "Sex", "Year")) %>%
  filter(mvmtClass != "resident") %>%
  mutate(indYr = paste0(ID, "_", Year)) %>%
  pull(indYr) %>%
  unique() %>%
  length()
ssfMigr = myDat %>%
  group_by(stepID) %>%
  filter(MvmtClass != "resident") %>%
  ungroup() %>%
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
MigrMovements_fullModel = plotSSF(ssfMigr) + 
  ylab(expression(paste("Model ", beta))) +
  geom_text(aes(x = "Elevation*Hot day (start)",
                y = 0.2),
            label = paste0("n = ", migrIndYrs),hjust=0,
            size = 2.5,
            col = "grey40") +
  ggtitle("Migrants")
MigrMovements_fullModel

## Residents
resiIndYrs = myDat %>%
  left_join(migPat %>% mutate(Year = springYr),
            by = c("ID", "Sex", "Year")) %>%
  filter(mvmtClass == "resident") %>%
  mutate(indYr = paste0(ID, "_", Year)) %>%
  pull(indYr) %>%
  unique() %>%
  length()
ssfResi = myDat %>%
  group_by(stepID) %>%
  filter(MvmtClass == "resident") %>%
  ungroup() %>%
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
ssfResi
ResiMovements_fullModel = plotSSF(ssfResi) + 
  ylab(expression(paste("Model ", beta))) +
  geom_text(aes(x = "Elevation*Hot day (start)",
                y = 0.2),
            label = paste0("n = ", resiIndYrs),hjust=0,
            size = 2.5,
            col = "grey40") +
  ggtitle("Residents")
ResiMovements_fullModel

## > Compile mig vs res plots ----
migrMvmtComparison = cowplot::plot_grid(
  MigrMovements_fullModel +
    theme(axis.text.x = element_blank()) + 
    theme(plot.background = element_rect(fill="transparent"),
          panel.background = element_rect(fill="transparent")),
  ResiMovements_fullModel +
    theme(plot.background = element_rect(fill="transparent"),
          panel.background = element_rect(fill="transparent")),
  rel_heights = c(0.5,1),
  align = "v",
  axis = "trbl",
  labels = c("A","B"),
  nrow = 2)
migrMvmtComparison

#### Variability by sex ----
## Ewes
eweIndYrs = myDat %>%
  left_join(migPat %>% mutate(Year = springYr),
            by = c("ID", "Sex", "Year")) %>%
  filter(Sex=="F") %>%
  mutate(indYr = paste0(ID, "_", Year)) %>%
  pull(indYr) %>%
  unique() %>%
  length()
ssfEwes = myDat %>%
  filter(Sex == "F") %>%
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
ewes_fullModel = plotSSF(ssfEwes) + 
  ylab(expression(paste("Model ", beta))) +
  geom_text(aes(x = "Elevation*Hot day (start)",
                y = 0.2),
            label = paste0("n = ", eweIndYrs),hjust=0,
            size = 2.5,
            col = "grey40") +
  ggtitle("Females")
ewes_fullModel

## Rams
ramIndYrs = myDat %>%
  left_join(migPat %>% mutate(Year = springYr),
            by = c("ID", "Sex", "Year")) %>%
  filter(Sex=="M") %>%
  mutate(indYr = paste0(ID, "_", Year)) %>%
  pull(indYr) %>%
  unique() %>%
  length()
ssfRams = myDat %>%
  filter(Sex == "M") %>%
  ungroup() %>%
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

ramPrecipSummary = as.data.frame(summary(ssfRams)$coefficients) %>%
  mutate("term" = rownames(.),
         log_coef = exp(coef),
         log_min = exp(coef - `se(coef)`),
         log_max = exp(coef + `se(coef)`)) %>%
  filter(term == "Precip_end_scaled")
ramPrecipEndpoints = c(ramPrecipSummary$log_min,
                       ramPrecipSummary$log_max)

rams_fullModel = plotSSF(ssfRams) + 
  geom_segment(aes(x = "Precipitation", xend = "Precipitation", 
                   y = ramPrecipEndpoints[1], yend = 2.95),
               col = "grey60") +
  geom_text(aes(x = "Precipitation", y = 3),
            label = round(ramPrecipEndpoints[2], 2),
            size = 2.5, col = "grey60", hjust = 1.2) +
  ylab(expression(paste("Model ", beta))) +
  geom_text(aes(x = "Elevation*Hot day (start)",
                y = 0.2),
            label = paste0("n = ", ramIndYrs),hjust=0,
            size = 2.5,
            col = "grey40") +
  ggtitle("Males")
rams_fullModel

cowplot::plot_grid(ewes_fullModel + theme(axis.text.x = element_blank()),
                   rams_fullModel + theme(axis.text.x = element_text(size = 6)),
                   nrow = 2,
                   rel_heights = c(0.6,1))


#### Variability by migration rate ----
directIndYrs = myDat %>%
  left_join(migPat %>% mutate(Year = springYr),
            by = c("ID", "Sex", "Year")) %>%
  filter(uphillRate <= 1) %>%
  mutate(indYr = paste0(ID, "_", Year)) %>%
  pull(indYr) %>%
  unique()%>%
  length()
ssfDirect = myDat %>%
  left_join(migPat %>% mutate(Year = springYr),
            by = c("ID", "Sex", "Year")) %>%
  filter(uphillRate <= 1) %>%
  ungroup() %>%
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
directMovers = plotSSF(ssfDirect) + 
  ylab(expression(paste("Model ", beta))) +
  ylim(0,3.5) +
  geom_text(aes(x = "Elevation*Hot day (start)",
                y = 0.2),
            label = paste0("n = ", fastIndYrs),hjust=0,
            size = 2.5,
            col = "grey40") +
  ggtitle("Direct migrants")
directMovers

slowIndYrs = myDat %>%
  left_join(migPat %>% mutate(Year = springYr),
            by = c("ID", "Sex", "Year")) %>%
  filter(uphillRate >= 7) %>%
  mutate(indYr = paste0(ID, "_", Year)) %>%
  pull(indYr) %>%
  unique() %>%
  length()
ssfSlow = myDat %>%
  left_join(migPat %>% mutate(Year = springYr),
            by = c("ID", "Sex", "Year")) %>%
  filter(uphillRate >= 7) %>%
  ungroup() %>%
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
slowMovers = plotSSF(ssfSlow) + 
  ylab(expression(paste("Model ", beta))) +
  ylim(0,3.5) +
  geom_text(aes(x = "Elevation*Hot day (start)",
                y = 0.2),
            label = paste0("n = ", slowIndYrs),hjust=0,
            size = 2.5,
            col = "grey40") +
  ggtitle("Slow migrants")
slowMovers

fastIndYrs = myDat %>%
  left_join(migPat %>% mutate(Year = springYr),
            by = c("ID", "Sex", "Year")) %>%
  filter(uphillRate < 7) %>%
  mutate(indYr = paste0(ID, "_", Year)) %>%
  pull(indYr) %>%
  unique() %>%
  length()
ssffast = myDat %>%
  left_join(migPat %>% mutate(Year = springYr),
            by = c("ID", "Sex", "Year")) %>%
  filter(uphillRate < 7) %>%
  ungroup() %>%
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
fastMovers = plotSSF(ssffast) + 
  ylab(expression(paste("Model ", beta))) +
  ylim(0,3.5) +
  geom_text(aes(x = "Elevation*Hot day (start)",
                y = 0.2),
            label = paste0("n = ", fastIndYrs),hjust=0,
            size = 2.5,
            col = "grey40") +
  ggtitle("Fast migrants")
fastMovers

cowplot::plot_grid(fastMovers + theme(axis.text.x = element_blank()),
                   slowMovers + theme(axis.text.x = element_text(size = 6)),
                   nrow = 2,
                   rel_heights = c(0.6,1))

## > Compile plots ----
myAxLim = ylim(0,3)
## Rams vs ewes
sexComparison = cowplot::plot_grid(
  ewes_fullModel +
    myAxLim +
    theme(axis.text.x = element_blank()) + 
    theme(plot.background = element_rect(fill="transparent"),
          panel.background = element_rect(fill="transparent")),
  rams_fullModel + 
    myAxLim +
    theme(axis.text.x = element_blank()) + 
    theme(plot.background = element_rect(fill="transparent"),
          panel.background = element_rect(fill="transparent")),
  align = "v",
  axis = "trbl",
  labels = c("A","B"),
  nrow = 1)
sexComparison
## Migrants vs residents
stratComparison = cowplot::plot_grid(
  ResiMovements_fullModel + 
    myAxLim +
    theme(axis.text.x = element_blank()) + 
    theme(plot.background = element_rect(fill="transparent"),
          panel.background = element_rect(fill="transparent")),
  MigrMovements_fullModel +
    myAxLim +
    theme(axis.text.x = element_blank()) + 
    theme(plot.background = element_rect(fill="transparent"),
          panel.background = element_rect(fill="transparent")),
  align = "v",
  axis = "trbl",
  labels = c("C","D"),
  nrow = 1)
stratComparison
## fast migrants vs slow migrants
speedComparison = cowplot::plot_grid(
  fastMovers +
    myAxLim +
    theme(plot.background = element_rect(fill="transparent"),
          panel.background = element_rect(fill="transparent")),
  slowMovers + 
    myAxLim +
    theme(plot.background = element_rect(fill="transparent"),
          panel.background = element_rect(fill="transparent")),
  align = "v",
  axis = "trbl",
  labels = c("E","F"),
  nrow = 1)
speedComparison


#### Combining migrants vs residents, and downhill vs uphill mvmts ----
subGroupMvmtComparison = cowplot::plot_grid(sexComparison,
                                            stratComparison,
                                            speedComparison,
                                            nrow = 3,
                                            rel_heights = c(1,1,2))

ggsave("plots/subgroupMvmtComparison.jpg",
       subGroupMvmtComparison,
       width = 6, height = 8, units = "in", dpi = 300)
ggsave("plots/subgroupMvmtComparison.png",
       subGroupMvmtComparison,
       width = 6, height = 8, units = "in", dpi = 300)

