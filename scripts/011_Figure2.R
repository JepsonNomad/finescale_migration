#### Load packages ----
library(tidyverse)
library(terra)
library(sf)
library(lme4)
library(lmerTest)

#### Define parameters ----
yearsOfInterest = c(2003:2022)

HUlookups = data.frame(FullNames = c("Olancha Peak", "Mt. Langley", 
                                     "Big Arroyo", "Laurel Creek", 
                                     "Mt. Williamson", "Mt.Baxter", 
                                     "Bubbs Creek", "Sawmill Canyon",
                                     "Taboose Creek", "Coyote Ridge",
                                     "Wheeler Ridge", "Convict Creek",
                                     "Mt. Gibbs", "Mt. Warren",
                                     "Green Creek", "Twin Lakes",
                                     "Cathedral Range", "Black Divide"),
                       CorrectNames = c("Olancha Peak", "Mt. Langley", 
                                        "Big Arroyo", "Laurel Creek", 
                                        "Mt. Williamson", "Mt. Baxter", 
                                        "Bubbs Creek", "Sawmill Canyon",
                                        "Taboose Creek", "Coyote Ridge",
                                        "Wheeler Ridge", "Convict Creek",
                                        "Mt. Gibbs", "Mt. Warren",
                                        "Green Creek", "Twin Lakes",
                                        "Cathedral Range", "Black Divide"),
                       Abbrevs = c("Ol","Ln","Ba","Lr","Wl","Bx","Bb","Sw",
                                   "Tb","Cy","Wh","Cn","Gb","Wr","Gc","Tw",
                                   "Ca","Bd"))

#### Import data ----
## > Migration phenology ----
migPhenol = read_csv("data/mig_PatternsCompiled.csv") %>%
  dplyr::rename(Individual = ID) %>%
  mutate(ID = paste0(Individual,"_",springYr)) %>%
  filter(springYr %in% yearsOfInterest)
## Tidy up the migration data
indParams = migPhenol %>%
  dplyr::select(ID, Sex, mvmtClass, upDOY, uphillRate, HerdUnit) %>%
  dplyr::rename(HU = HerdUnit)
indParams$HerdUnit = HUlookups$FullNames[match(indParams$HU,
                                               HUlookups$Abbrevs)]
indParams = indParams %>%
  mutate(HerdUnit = case_when(HerdUnit == "Mt.Baxter" ~ "Mt. Baxter",
                              TRUE ~ HerdUnit))

## > Landscape phenology ----
## Find raster data
snowMeltFPs = list.files(path = "../SSN/snowPhenologyMosaics",
                         pattern  ="*.tif", full.names = TRUE)
greenUpFPs = list.files(path = "../MOD13Q1/data/QA1",
                        pattern = "phenology", full.names = TRUE)

## Note: snowfall date is layer 1, snowmelt date is layer 2
snowMelt = lapply(snowMeltFPs,
                  function(x){
                    r = rast(x,lyr=2)
                    r[r>900] <- NA
                    r = r - 227
                    return(r)
                  }) %>%
  rast()
names(snowMelt) = sapply(str_split(sapply(str_split(snowMeltFPs, "/"),
                                          "[[", 5), ".tif"),
                         "[[", 1)

## Note: greenup date is layer 1
greenUp = lapply(greenUpFPs,
                 function(x){
                   rast(x,lyr=1)
                 }) %>%
  rast()
names(greenUp) = sapply(str_split(sapply(str_split(greenUpFPs, "_"),
                                         "[[", 2), ".tif"),
                        "[[", 1)

## > Elevation ----
DEM = rast("../USGS_3DEP/data/DEM_UTM.tif") %>%
  resample(greenUp)

## > Herd units  ----
HUs = st_read("../SNBS_HerdUnits/Cathedral/HerdUnits_UpdatedMetadata.shp") %>%
  st_transform(st_crs(snowMelt[[1]])) %>%
  mutate(NAME = case_when(NAME == "Mt.Baxter" ~ "Mt. Baxter",
                          TRUE ~ NAME))

#### Data wrangling ----
## Make a resampled snowmelt layer to match greenUp:
snowMeltResampled = snowMelt %>%
  resample(greenUp)

## Filter greenup and snowmelt to years of interest
greenUpTimeframe = greenUp[[names(greenUp) %in% yearsOfInterest]]
snowMeltTimeframe = snowMeltResampled[[names(snowMeltResampled) %in% yearsOfInterest]]


## Extract snow and greenup phenology data by herd unit
## Pull out mean and standard deviation for each herd unit
snExt = extract(snowMeltTimeframe, vect(HUs), fun = "mean", na.rm = T)
guExt = extract(greenUpTimeframe, vect(HUs), fun = "mean", na.rm = T)
snSD = extract(snowMeltTimeframe, vect(HUs), fun = "sd", na.rm = T)
guSD = extract(greenUpTimeframe, vect(HUs), fun = "sd", na.rm = T)

## Add a herd unit column to phenology dataesets
snExt$HerdUnit = HUs$NAME
guExt$HerdUnit = HUs$NAME
snSD$HerdUnit = HUs$NAME
guSD$HerdUnit = HUs$NAME

## Pivot landscape phenology datasets to long format for joining
snLong = snExt %>%
  dplyr::select(-ID) %>%
  pivot_longer(-HerdUnit,
               values_to = "Snowmelt") %>%
  mutate(Year = as.numeric(gsub("X","",name))) %>%
  dplyr::select(HerdUnit, Year, Snowmelt)
guLong = guExt %>%
  dplyr::select(-ID) %>%
  pivot_longer(-HerdUnit,
               values_to = "Greenup") %>%
  mutate(Year = as.numeric(gsub("X","",name))) %>%
  dplyr::select(HerdUnit, Year, Greenup)
snSDLong = snSD %>%
  dplyr::select(-ID) %>%
  pivot_longer(-HerdUnit,
               values_to = "SnowmeltSD") %>%
  mutate(Year = as.numeric(gsub("X","",name))) %>%
  dplyr::select(HerdUnit, Year, SnowmeltSD)
guSDLong = guSD %>%
  dplyr::select(-ID) %>%
  pivot_longer(-HerdUnit,
               values_to = "GreenupSD") %>%
  mutate(Year = as.numeric(gsub("X","",name))) %>%
  dplyr::select(HerdUnit, Year, GreenupSD)

## Compile all landscape phenology data
mydata = full_join(snLong,
                   guLong,
                   by = c("HerdUnit",
                          "Year")) %>%
  left_join(snSDLong,
            by = c("HerdUnit",
                   "Year")) %>%
  left_join(guSDLong,
            by = c("HerdUnit",
                   "Year"))
## Preview
mydata %>%
  ggplot() +
  geom_line(aes(x = Year,
                y = Snowmelt,
                group = HerdUnit),
            col = "lightblue") +
  geom_line(aes(x = Year,
                y = Greenup,
                group = HerdUnit),
            col = "green3") +
  facet_wrap(~HerdUnit) +
  ylab("Day of year") +
  CJsBasics::BasicTheme

#### Migration phenology in response to environmental conditions
migdat = indParams %>%
  separate(ID, into = c("IndID","Year"), sep = "_") %>%
  mutate(Year = as.numeric(Year)) %>%
  left_join(mydata, by = c("Year","HerdUnit")) %>% 
  filter(!is.na(mvmtClass))

migdat %>%
  arrange(Year)

## > Overall landscape phenology paragraph ----
# Green-up and snowmelt timing varied from year to year and both advanced significantly over the course of the study (Figure2A-C). A linear mixed effects model with green-up or snowmelt timing as the response variable, year as the predictor, and herd unit as a random effect, revealed that snowmelt timing advanced at a rate of 8.03±1.80 days per decade, and green-up timing advanced at a rate of 6.14±1.57 days per decade (p<0.001 in both cases). Variability within herd units (indexed as the standard deviation of green-up or snowmelt timing each year) did not significantly change over the course of the time series (p > 0.05). 
ts1Data = mydata %>%
  mutate(HerdUnit = factor(
    HerdUnit, 
    levels = HUlookups$CorrectNames[nrow(HUlookups):1])) %>%
  filter(HerdUnit %in% migdat$HerdUnit) %>%
  filter(Year %in% yearsOfInterest)
ma = lmer(Greenup ~ Year +
            (1 | HerdUnit),
          data = ts1Data)
summary(ma)
mb = lmer(Snowmelt ~ Year +
            (1 | HerdUnit),
          data = ts1Data)
summary(mb)
ma2 = lmer(GreenupSD ~ Year +
             (1 | HerdUnit),
           data = ts1Data)
summary(ma2)
mb2 = lmer(SnowmeltSD ~ Year +
             (1 | HerdUnit),
           data = ts1Data)
summary(mb2)

ts1 = ts1Data %>%
  ggplot() +
  annotate(x = 2003.5, y = 250, label = "Green-up",
           geom = "text",
           hjust = 0) +
  annotate(x = 2003.5, y = -70, label = "Snowmelt",
           geom = "text",
           hjust = 0) +
  geom_ribbon(aes(x = Year,
                  ymin = Greenup - GreenupSD,
                  ymax = Greenup + GreenupSD,
                  fill = HerdUnit,
                  group = HerdUnit),
              alpha = 0.1,
              show.legend = F) +
  geom_ribbon(aes(x = Year,
                  ymin = Snowmelt - SnowmeltSD,
                  ymax = Snowmelt + SnowmeltSD,
                  fill = HerdUnit,
                  group = HerdUnit),
              alpha = 0.1,
              show.legend = F) +
  geom_line(aes(x = Year,
                y = Greenup,
                col = HerdUnit,
                group = HerdUnit),
            lwd = 0.35) +
  geom_line(aes(x = Year,
                y = Snowmelt,
                col = HerdUnit,
                group = HerdUnit),
            lty = 2,
            lwd = 0.35) +
  scale_color_manual("Herd unit",
                     values = CJsBasics::KellyCols[c(2:20)]) +
  scale_fill_manual("Herd unit",
                    values = CJsBasics::KellyCols[c(2:20)]) +
  xlab("Year") +
  ylab("Day of year") +
  CJsBasics::BasicTheme

ts2 = ts1Data %>%
  group_by(Year) %>%
  summarize(Greenup = mean(Greenup),
            Snowmelt = mean(Snowmelt)) %>%
  ggplot() +
  geom_point(aes(x = Year, y = Snowmelt)) +
  stat_smooth(aes(x = Year, y = Snowmelt),
              method = "lm",
              col = "grey60",
              linetype = 2,
              se = F) +
  ggtitle("Snowmelt timing") +
  ylab("Day of year") + 
  CJsBasics::BasicTheme
ts3 = ts1Data %>%
  group_by(Year) %>%
  summarize(Greenup = mean(Greenup)) %>%
  ggplot() +
  geom_point(aes(x = Year, y = Greenup)) +
  stat_smooth(aes(x = Year, y = Greenup),
              method = "lm",
              col = "grey60",
              se = F) +
  ggtitle("Green-up timing") +
  ylab("Day of year") +
  CJsBasics::BasicTheme

## Compile figures 2a-c
ts1bc = cowplot::plot_grid(ts2 + scale_x_continuous(breaks=c(2005,2020)) +
                             theme(plot.margin = margin(0.1,0.25,0.1,0,"in")), 
                           ts3 + scale_x_continuous(breaks=c(2005,2020)) +
                             theme(plot.margin = margin(0.1,0.25,0.1,0,"in")),
                           ncol = 2,
                           align = "hv",
                           labels = c("B","C"))
ts1abc = cowplot::plot_grid(ts1 +
                              theme(legend.key.height = unit(0.125, "in"),
                                    legend.spacing.y = unit(0.025, "in"),
                                    legend.background = element_blank()),
                            ts1bc, 
                            ncol = 2,
                            rel_widths = c(1,0.8),
                            labels = c("A",""))
## > Green-up timing paragraph ----
# In a mixed-effects model with green-up timing as the response variable, snowmelt timing as the predictor variable, and year as a random intercept, green-up was 7.7±0.15 days later per 10-day delay in snowmelt (p < 0.001; conditional R2 = 0.92)
m1 = lmer(Greenup ~ Snowmelt +
            (1 | HerdUnit),
          data = migdat)
summary(m1)
# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: Greenup ~ Snowmelt + (1 | HerdUnit)
#    Data: migdat
# 
# REML criterion at convergence: 4769.6
# 
# Scaled residuals: 
#      Min       1Q   Median       3Q      Max 
# -2.50466 -0.71451 -0.02864  0.84624  2.32874 
# 
# Random effects:
#  Groups   Name        Variance Std.Dev.
#  HerdUnit (Intercept) 46.80    6.841   
#  Residual             48.49    6.963   
# Number of obs: 702, groups:  HerdUnit, 14
# 
# Fixed effects:
#              Estimate Std. Error        df t value Pr(>|t|)    
# (Intercept) 141.68689    1.88869  13.83904   75.02   <2e-16 ***
# Snowmelt      0.77372    0.01462 692.39508   52.91   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#          (Intr)
# Snowmelt -0.171
MuMIn::r.squaredGLMM(m1)
# R2m      R2c
# [1,] 0.834151 0.915611

# The mean difference in timing between peak snowmelt and peak green-up was 136.7±0.86 days
grnsnowTtest = t.test(migdat$Greenup, migdat$Snowmelt, paired = T)
grnsnowTtest$estimate
# mean of the differences 
# 136.7478 
grnsnowTtest$estimate - grnsnowTtest$conf.int[1]
# 0.8643924 


## .. Figure 2c ----
## Create regression lines from mixed effects model
m1predData = mydata %>%
  filter(HerdUnit %in% migdat$HerdUnit) %>%
  arrange(Snowmelt)
m1Preds = predict(m1, newdata = m1predData, interval = "prediction")
## Refactor the data and add in the regression lines; select relevant years
p1data = mydata %>%
  mutate(HerdUnit = factor(
    HerdUnit, 
    levels = HUlookups$CorrectNames[nrow(HUlookups):1])) %>%
  filter(HerdUnit %in% migdat$HerdUnit) %>%
  arrange(Snowmelt) %>%
  mutate(colm1 = m1Preds) %>%
  filter(Year %in% yearsOfInterest)
## Plot data
p1 = p1data %>%
  ggplot(aes(x = Snowmelt, y = Greenup)) +
  geom_point(aes(col = HerdUnit),
             alpha = 0.75) +
  geom_line(aes(y = colm1,
                col = HerdUnit)) +
  # geom_abline() +
  scale_color_manual("Herd unit",
                     values = CJsBasics::KellyCols[c(2:20)]) +
  xlab("Snowmelt timing (day of year)") +
  ylab("Green-up timing (day of year)") +
  coord_equal() +
  CJsBasics::BasicTheme +
  theme(legend.key.height = unit(0.125, "in"),
        legend.spacing.y = unit(0.025, "in"),
        legend.background = element_blank(),
        legend.position = c(0.825,0.305))

## > Migration timing paragraph ----
# Mean migration timing occurred before mean green-up timing at the herd unit level in 75.2% of cases, however anomalously late migrations were observed in several cases when bighorn moved uphill as late as mid-August 
migdat %>%
  group_by(HerdUnit, Year, Greenup) %>%
  summarize(upDOY = mean(upDOY, na.rm = T)) %>%
  ungroup() %>%
  mutate(upDiff = upDOY < Greenup) %>%
  filter(!is.nan(upDOY)) %>%
  pull(upDiff) %>%
  sum(.)
## 85
migdat %>%
  group_by(HerdUnit, Year, Greenup) %>%
  summarize(upDOY = mean(upDOY, na.rm = T)) %>%
  ungroup() %>%
  mutate(upDiff = upDOY < Greenup) %>%
  filter(!is.nan(upDOY)) %>%
  pull(upDiff) %>%
  length()
## 113
85/113
# [1] 0.7522124

migdat %>%
  group_by(HerdUnit, Year, Greenup) %>%
  summarize(upDOY = mean(upDOY, na.rm = T)) %>%
  mutate(upDiff = upDOY > Greenup) %>%
  filter(upDiff) %>%
  arrange(-(upDOY))
## DOY 228 = Aug 16

# The earliest migration (relative to green-up timing) occurred 124 days prior to mean green-up timing (at the herd unit level), while the latest was 79 days after mean green-up timing
migdat %>% 
  group_by(HerdUnit, Year, Greenup) %>%
  summarize(upDOY = mean(upDOY, na.rm = T)) %>%
  mutate(upDiff = upDOY - Greenup) %>% 
  pull(upDiff) %>%
  range(na.rm=T) %>%
  round()

# A linear mixed-effects model with migration timing as the response variable, green-up timing as a predictor, and herd unit identity and year as random intercepts, revealed that the midpoint of migration was 4.0±1.1 days later per 10-day delay in green-up timing (p = 0.002; conditional R2 = 0.46)
m2 = lmer(upDOY ~ Greenup + Sex +
            (1 | HerdUnit) + (1 | Year),
          data = migdat)
summary(m2)

# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: upDOY ~ Greenup + Sex + (1 | HerdUnit) + (1 | Year)
#    Data: migdat
# 
# REML criterion at convergence: 3621.5
# 
# Scaled residuals: 
#     Min      1Q  Median      3Q     Max 
# -4.2663 -0.4918 -0.0391  0.4641  4.6223 
# 
# Random effects:
#  Groups   Name        Variance Std.Dev.
#  Year     (Intercept)  15.49    3.936  
#  HerdUnit (Intercept) 535.53   23.142  
#  Residual             760.52   27.578  
# Number of obs: 379, groups:  Year, 20; HerdUnit, 14
# 
# Fixed effects:
#             Estimate Std. Error      df t value Pr(>|t|)    
# (Intercept)   85.131     18.012  19.148   4.726 0.000144 ***
# Greenup        0.406      0.105  16.483   3.866 0.001304 ** 
# SexM           5.915      3.471 365.090   1.704 0.089223 .  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#         (Intr) Greenp
# Greenup -0.928       
# SexM    -0.104  0.054
MuMIn::r.squaredGLMM(m2)
# R2m       R2c
# [1,] 0.07037181 0.4609417

## Difference in estimated random intercept between rams and ewes:
sexDiffs = ranef(m2)$Sex
sexDiffs[rownames(sexDiffs)=="F",] - sexDiffs[rownames(sexDiffs)=="M",]
ggplot() +
  geom_histogram(aes(x = ranef(m2)$Year[,1]))
ggplot() +
  geom_histogram(aes(x = ranef(m2)$HerdUnit[,1]))

## .. Figure 2b ----
p2pltdata = migdat %>%
  group_by(HerdUnit, Year, Greenup) %>%
  summarize(count = sum(!is.na(upDOY)),
            variability = sd(upDOY, na.rm = T),
            se = variability/sqrt(count),
            upDOY = mean(upDOY, na.rm = T)) %>%
  ungroup() %>%
  filter(count > 2) %>%
  group_by(HerdUnit) %>%
  mutate(HUYrCounts = sum(!is.na(count))) %>%
  ungroup() %>%
  mutate(HerdUnit = factor(
    HerdUnit,
    levels = HUlookups$CorrectNames[nrow(HUlookups):1])) %>%
  filter(HUYrCounts > 4)
## Get colors sorted with first panel
p2cols = match(unique(sort(p2pltdata$HerdUnit)), 
               sort(unique(p1data$HerdUnit)))

p2 = p2pltdata %>%
  ggplot() +
  geom_abline() +
  geom_errorbar(aes(x = Greenup, 
                    ymin = c(upDOY-se),
                    ymax = c(upDOY+se),
                    col = HerdUnit),
                alpha = 0.5) +
  geom_point(aes(x = Greenup, y = upDOY,
                 size = count),
             col = "white") +
  geom_point(aes(x = Greenup, y = upDOY,
                 size = count,
                 col = HerdUnit),
             alpha = 0.75) + 
  stat_smooth(aes(x = Greenup, y = upDOY,
                  group = HerdUnit,
                  col = HerdUnit),
              method = "lm",
              fill = "#e1e1e190",
              se = F,
              show.legend = F) +
  xlab("Green-up timing (day of year)") +
  ylab("Migration timing (day of year)") +
  scale_size_continuous("Count", 
                        breaks = c(3,6,12)) +
  scale_alpha_continuous("Variance", 
                         trans = "reverse") +
  scale_color_manual("Herd unit",
                     values = CJsBasics::KellyCols[c(1+p2cols)]) +
  coord_equal() +
  xlim(c(50,
         max(c(p2pltdata$upDOY+p2pltdata$se, 
               p2pltdata$Greenup, 212),
             na.rm = T))) +
  ylim(c(50,
         max(c(p2pltdata$upDOY+p2pltdata$se,
               p2pltdata$Greenup, 212),
             na.rm = T))) +
  CJsBasics::BasicTheme +
  theme(legend.key.height = unit(0.125, "in"),
        legend.background = element_blank(),
        legend.margin = margin(),
        legend.spacing.y = unit(0.025, "in"),
        legend.position = c(0.2,0.775),
        legend.title.align = 0) +
  guides(color = guide_legend(order = 1),
         size = guide_legend(order= 2))
p2

## .. Compile Figures 2d and 2e ----
## Figure 2
p3 = cowplot::plot_grid(p1,p2,
                        nrow = 1,
                        labels = c("D","E"),
                        align = "h")
p3

fig2 = cowplot::plot_grid(ts1abc,
                          p3,
                          rel_heights = c(0.6,1),
                          axis = "trbl",
                          nrow = 2)
ggsave("plots/figure2.jpg",
       fig2,
       width = 9, height = 7, units = "in", dpi = 300)
