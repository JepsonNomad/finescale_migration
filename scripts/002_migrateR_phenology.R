startTime = Sys.time()
#### Extracting migration phenology from prepped data
#### Christian John
#### 29 Nov 2023
## Takes about 5 sec to run

#### Load packages ----
library(plyr)
library(tidyverse)
library(lme4)
library(lmerTest)

source("scripts/000_samplingParameters.R")

#### Import data ----
dat = read_csv("data/migrationPatterns_allInds.csv")

#### Data wrangling ----
# Create empty columns for date calculations
dat$downhillMvmt = as.Date(NA)
dat$uphillMvmt = as.Date(NA)
# Create columns for range and migration parameters 
dat = dat %>%
  filter(springYr %in% myYears,
         !is.na(mvmtClass)) %>%
  mutate(
    highElev = case_when(delta < 0 ~ 
                           startElev,
                         TRUE ~ endElev),
    lowElev = case_when(delta > 0 ~ 
                          startElev,
                        TRUE ~ endElev),
    downhillMvmt = case_when(delta < 0 ~ 
                               migTiming1,
                             TRUE ~ migTiming2),
    uphillMvmt = case_when(delta > 0 ~ 
                             migTiming1,
                           TRUE ~ migTiming2),
    downhillRate = case_when(delta < 0 ~ 
                               phi,
                             TRUE ~ phi2),
    uphillRate = case_when(delta > 0 ~ 
                             phi,
                           TRUE ~ phi2),
    downDOY = as.numeric(strftime(downhillMvmt, 
                                  format = "%j")),
    upDOY = as.numeric(strftime(uphillMvmt,
                                format = "%j")))

#### Summaries ----
dat %>%
  group_by(springYr, HerdUnit) %>%
  summarize(nInds = n()) %>%
  print(n = 300)

## How many individuals
dat %>%
  select(ID) %>%
  unique()
# A tibble: 311 × 1

## How many animal-years
dat %>%
  nrow() 
# 702

## How many inds per year
dat %>%
  group_by(springYr) %>%
  summarize(nInds = n()) %>%
  ungroup() %>%
  summarize(minInds = min(nInds),
            maxInds = max(nInds),
            meanInds = mean(nInds),
            seInds = plotrix::std.error(nInds))
# minInds maxInds meanInds seInds
# <int>   <int>    <dbl>  <dbl>
#   1       1      74     35.1   5.13

## How many herd units
dat %>%
  select(HerdUnit) %>%
  unique() %>%
  nrow()
# [1] 14

## "Supplementary materials S3"
dat %>%
  filter(!(mvmtClass != "resident" & is.na(uphillRate))) %>%
  mutate(Speed = ifelse(uphillRate < 7, "fast","slow"),
         Strategy = ifelse(mvmtClass == "resident",
                           "resident",
                           "migrant"),
         Year = springYr) %>%
  mutate(Speed = replace_na(Speed, replace = "")) %>%
  group_by(Year, HerdUnit, Strategy, Speed) %>%
  summarize(nInds = n()) %>%
  ungroup() %>%
  mutate(nInds = as.character(nInds)) %>%
  pivot_wider(names_from = "Year", values_from = "nInds",
              values_fill = "") %>%
  arrange(HerdUnit, Strategy, Speed) %>%
  rename("Herd unit" = HerdUnit) %>%
  kableExtra::kable(table.attr = "style='width:100%;'") %>%
  kableExtra::save_kable(file = "plots/mvmtClass_breakdown_byYear.jpg")

## "Accordingly, we observed a high rate of switching between migration and residency by individuals among years (mean strategy switch rate = 0.46±0.04) and switching between fast and slow migration tactics among years (mean rate switch rate = 0.62±0.06)"
classRepeatability = dat %>%
  arrange(springYr) %>%
  group_by(ID) %>%
  mutate(indYRs = n()) %>%
  ungroup() %>%
  filter(indYRs > 1) %>%
  group_by(ID) %>%
  mutate(classMatch = mvmtClass == lag(mvmtClass)) %>%
  ungroup() %>%
  arrange(ID, springYr) %>%
  select(ID, springYr, classMatch)
classRepeatability %>%
  filter(!is.na(classMatch)) %>%
  group_by(springYr) %>%
  summarize(retention = sum(classMatch, na.rm = T) / n()) %>%
  summarise(mean(retention))
classRepeatability %>%
  filter(!is.na(classMatch)) %>%
  group_by(springYr) %>%
  summarize(retention = sum(classMatch, na.rm = T) / n(),
                              digits = 3) %>%
  summarise(plotrix::std.error(retention))

speedRepeatability = dat %>%
  filter(!is.na(uphillRate)) %>%
  mutate(mvmtspeed = ifelse(uphillRate >= 7,
                            "slow",
                            "fast")) %>%
  arrange(springYr) %>%
  group_by(ID) %>%
  mutate(indYRs = n()) %>%
  ungroup() %>%
  filter(indYRs > 1) %>%
  group_by(ID) %>%
  mutate(speedMatch = mvmtspeed == lag(mvmtspeed)) %>%
  ungroup() %>%
  arrange(ID, springYr) %>%
  dplyr::select(ID, springYr, speedMatch)
speedRepeatability %>%
  filter(!is.na(speedMatch)) %>%
  group_by(springYr) %>%
  summarize(retention = sum(speedMatch, na.rm = T) / n()) %>%
  summarise(mean(retention))
speedRepeatability %>%
  filter(!is.na(speedMatch)) %>%
  group_by(springYr) %>%
  summarize(retention = sum(speedMatch, na.rm = T) / n(),
            digits = 3) %>%
  summarise(plotrix::std.error(retention))


#### Save output ----
write_csv(dat, 
          file = "data/mig_PatternsCompiled.csv")

#### Epilogue ----
stopTime = Sys.time()
print(stopTime - startTime)
sessionInfo()
