library(sf)
library(tidyverse)
library(terra)

#### Figure 1 ----
indsOfInterest =  c("S107", "S55", "S89")

## Migration summary data
migPhenol = read_csv("data/mig_PatternsCompiled.csv") %>%
  rename(Individual = ID) %>%
  mutate(ID = paste0(Individual,"_",springYr))
indParams = migPhenol %>%
  select(ID, Sex, mvmtClass, upDOY, uphillRate) %>%
  mutate(mvmtClass = ifelse(mvmtClass == "resident",
                            "resident",
                            "migrant"))

## Raster data
DEM = rast("../USGS_3DEP/data/DEM_UTM.tif")

## Vector data
sheep_sf = st_read(
  "../data_ts/SNBS_Locs_NoTranslocated/AllSierra_NoTranslocatedfor1yr.shp") %>%
  rename(ID = AnimalID)
mysheep = sheep_sf %>%
  filter(ID %in% indsOfInterest,
         Year == 2008) %>%
  mutate(datetime = lubridate::ymd_hms(
    paste0(Date," ",Time))) %>%
  filter(as.Date(datetime) > as.Date("2008-02-15") &
           as.Date(datetime) < as.Date("2008-08-15"))
mysheep_vect = vect(mysheep)

#### Extract data ----
mysheep$DEM = terra::extract(DEM, mysheep_vect)[,2]

#### View modeled elevation ----
migPredictor = function(indID,
                        year = 2008){
  dfLine = migPhenol %>%
    filter(Individual == indID,
           springYr == year)
  
  burstYr = dfLine$burstYr
  mvmtcl = dfLine$mvmtClass
  if(is.na(mvmtcl)){
    outDF = data.frame("Date" = as.Date(NA),
                       "burstYr" = burstYr,
                       "predicted" = NA)
  }else{
    startDate = dfLine$indYrStartDate
    endDate = dfLine$indYrLastObs
    nDays = difftime(endDate,
                     startDate,
                     units = "days")
    dday = seq(from=1, to=as.numeric(nDays), by=1)
    
    predDates = seq.Date(from=as.Date(startDate),
                         to=as.Date(endDate - lubridate::days(1)),
                         by = "1 day")
    
    if(mvmtcl=="migrant"){
      delta = dfLine$delta
      theta = dfLine$theta
      phi = dfLine$phi
      rho = dfLine$rho
      theta2 = dfLine$theta2
      phi2 = dfLine$phi2
      gam = dfLine$gamma
      y = (delta/(1+exp((theta-dday)/phi))) +
        (-delta/(1+exp((theta+2*phi+2*phi2+rho-dday)/phi2))) + gam
    }else if(mvmtcl=="disperser"){
      delta = dfLine$delta
      theta = dfLine$theta
      phi = dfLine$phi
      gam = dfLine$gamma
      y = (delta/(1+exp((theta-dday)/phi))) + gam
    }else{
      gam = dfLine$gamma
      y = rep(gam, length(predDates))
    }
    outDF = data.frame("Date" = predDates,
                       "burstYr" = burstYr,
                       "predicted" = y,
                       "IndID" = indID)
  }
  return(outDF)
}

elevPredictions = do.call("rbind",
                          lapply(indsOfInterest,
                                 function(x){
                                   migPredictor(x)
                                 })) %>%
  filter(Date > as.Date("2008-02-15") &
           Date < as.Date("2008-08-15"))


#### Visualize ----
elevProfile = mysheep %>%
  st_drop_geometry() %>%
  ungroup() %>%
  ggplot() +
  geom_line(aes(x = datetime,
                y = DEM,
                col = ID),
            lwd = 0.5) +
  geom_line(data = elevPredictions,
            aes(x = as.POSIXct(Date),
                y = predicted,
                col = IndID),
            lwd = 2,
            alpha = 0.5) +
  xlab("Date") +
  ylab("Elevation (meters)") +
  scale_color_manual(values = wesanderson::wes_palette("Darjeeling2")[c(2,3,5)]) +
  CJsBasics::BasicTheme +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "transparent", 
                                        colour = NA),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_rect(fill='transparent')
  )

ggsave(x = elevProfile,
       file = "plots/elevProfileExample.png",
       bg = "transparent",
       width = 3, height = 3, units = "in", dpi = 300)
pdf("plots/elevProfileExample.pdf",
    width = 4, height = 3)
print(elevProfile)
dev.off()


