# Name: 2-freshwater-consumption-and-withdrawal.R 
# Description: Calculates total freshwater consumption for the year 2010 from Huang et al. (2018) across all sectors, reconciling the various agricultural estimates by using the median value of all alternatives.

# Use 'here' package for easy path management.  
library(here)

# Import all setup and user-defined functions in R/setup and R/udfs folders
invisible(sapply(paste0(here("R/setup"), "/", list.files(here("R/setup"))), source)) 
invisible(sapply(paste0(here("R/udf"), "/", list.files(here("R/udf"))), source))

# NOTE: While only the process for deriving the consumption data is shown below, the process is repeated for withdrawal data with no other changes made to the script.

setwd("D:/!! Geodatabase/Huang_WaterUse/2d netcdf")

# Decompress all data from 2d netcdf to 3d netcdf

# Domestic
Decompress2dto3dNCDf(Raw2dNCDF = ncdf4::nc_open("./cons_dom.nc"), 
                     varname = "cons_dom", 
                     var.unit = "mm/mo", 
                     CRS.set = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs", 
                     WriteLocation = here::here("Data", "WaterUse", "Consumption", "Domestic.nc"))

# Electric
Decompress2dto3dNCDf(Raw2dNCDF = ncdf4::nc_open("./cons_elec.nc"), 
                     varname = "cons_elec", 
                     var.unit = "mm/mo", 
                     CRS.set = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs", 
                     WriteLocation = here::here("Data", "WaterUse", "Consumption", "Electric.nc"))

# Irrigation H08
Decompress2dto3dNCDf(Raw2dNCDF = ncdf4::nc_open("./cons_irr_h08.nc"), 
                     varname = "cons_irr", 
                     var.unit = "mm/mo", 
                     CRS.set = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs", 
                     WriteLocation = here::here("Data", "WaterUse", "Consumption", "Irr_cons_h08.nc"))


# Irrigation lpjml
Decompress2dto3dNCDf(Raw2dNCDF = ncdf4::nc_open("./cons_irr_jpjml.nc"), 
                     varname = "cons_irr", 
                     var.unit = "mm/mo", 
                     CRS.set = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs", 
                     WriteLocation = here::here("Data", "WaterUse", "Consumption", "Irr_cons_lpjml.nc"))

# Irrigation pcr
Decompress2dto3dNCDf(Raw2dNCDF = ncdf4::nc_open("./cons_irr_pcr.nc"), 
                     varname = "cons_irr", 
                     var.unit = "mm/mo", 
                     CRS.set = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs", 
                     WriteLocation = here::here("Data", "WaterUse", "Consumption", "Irr_cons_pcr.nc"))

# Irrigation watergap
Decompress2dto3dNCDf(Raw2dNCDF = ncdf4::nc_open("./cons_irr_watergap.nc"), 
                     varname = "cons_irr", 
                     var.unit = "mm/mo", 
                     CRS.set = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs", 
                     WriteLocation = here::here("Data", "WaterUse", "Consumption", "Irr_cons_wgap.nc"))

# Livestock
Decompress2dto3dNCDf(Raw2dNCDF = ncdf4::nc_open("./cons_liv.nc"), 
                     varname = "cons_liv", 
                     var.unit = "mm/mo", 
                     CRS.set = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs", 
                     WriteLocation = here::here("Data", "WaterUse", "Consumption", "Livestock.nc"))

# Manufacturing
Decompress2dto3dNCDf(Raw2dNCDF = ncdf4::nc_open("./cons_mfg.nc"), 
                     varname = "cons_mfg", 
                     var.unit = "mm/mo", 
                     CRS.set = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs", 
                     WriteLocation = here::here("Data", "WaterUse", "Consumption", "Manufacturing.nc"))

# Mining
Decompress2dto3dNCDf(Raw2dNCDF = ncdf4::nc_open("./cons_min.nc"), 
                     varname = "cons_min", 
                     var.unit = "mm/mo", 
                     CRS.set = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs", 
                     WriteLocation = here::here("Data", "WaterUse", "Consumption", "Mining.nc"))


# Loop through each sector and sum over year 2010

# Initiate raster for each sector
domestic <- raster(ext = extent(-180, 180, -90, 90),
                   res=  c(0.5, 0.5),
                   crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
electric <- raster(domestic)
irr_H08 <- raster(domestic)
irr_lpjml <- raster(domestic)
irr_pcr <- raster(domestic)
irr_wgap <- raster(domestic)
livestock <- raster(domestic)
manufac <- raster(domestic)
mining <- raster(domestic)

domestic[] <- 0
electric[] <- 0
irr_H08[] <- 0
irr_lpjml[] <- 0
irr_pcr[] <- 0
irr_wgap[] <- 0
livestock[] <- 0
manufac[] <- 0
mining[] <- 0

for (i in 1:12) {
  Dom_m <- raster::stack(here("Data", "WaterUse", "Consumption", "Domestic.nc"))[[468+i]] 
  Elec_m <- raster::stack(here("Data", "WaterUse", "Consumption", "Electric.nc"))[[468+i]]
  Irr_h08_m <- raster::stack(here("Data", "WaterUse", "Consumption", "Irr_cons_h08.nc"))[[468+i]]
  Irr_lpj_m <- raster::stack(here("Data", "WaterUse", "Consumption", "Irr_cons_lpjml.nc"))[[468+i]]
  Irr_pcr_m <- raster::stack(here("Data", "WaterUse", "Consumption", "Irr_cons_pcr.nc"))[[468+i]]
  Irr_wgap_m <- raster::stack(here("Data", "WaterUse", "Consumption", "Irr_cons_wgap.nc"))[[468+i]]
  Liv_m <- raster::stack(here("Data", "WaterUse", "Consumption", "Livestock.nc"))[[468+i]]
  Mfg_m <- raster::stack(here("Data", "WaterUse", "Consumption", "Manufacturing.nc"))[[468+i]]
  Min_m <- raster::stack(here("Data", "WaterUse", "Consumption", "Mining.nc"))[[468+i]]
  
  Dom_m[is.na(Dom_m) | Dom_m < 0] <- 0
  Elec_m[is.na(Elec_m) | Elec_m < 0] <- 0
  Irr_h08_m[is.na(Irr_h08_m) | Irr_h08_m < 0] <- 0
  Irr_lpj_m[is.na(Irr_lpj_m) | Irr_lpj_m < 0] <- 0
  Irr_pcr_m[is.na(Irr_pcr_m) | Irr_pcr_m < 0] <- 0
  Irr_wgap_m[is.na(Irr_wgap_m) | irr_wgap < 0] <- 0
  Liv_m[is.na(Liv_m) | Liv_m < 0] <- 0
  Mfg_m[is.na(Mfg_m) | Mfg_m < 0] <- 0
  Min_m[is.na(Min_m) | Min_m < 0] <- 0
  
  domestic <- domestic + Dom_m
  electric <- electric + Elec_m
  irr_H08 <- irr_H08 + Irr_h08_m
  irr_lpjml <- irr_lpjml + Irr_lpj_m
  irr_pcr <- irr_pcr + Irr_pcr_m
  irr_wgap <- irr_wgap + Irr_wgap_m
  livestock <- livestock + Liv_m
  manufac <- manufac + Mfg_m
  mining <- mining + Min_m
  
  print(468+i)
}

# Calculate median irrigation consumption
Irr_med <- raster::stack(irr_H08, irr_lpjml, irr_pcr, irr_wgap)
Irr_med <- calc(Irr_med, median)

# Sum across all sectors
Total_cons <- raster::stack(domestic, electric, Irr_med, livestock, manufac, mining)
Total_cons <- calc(Total_cons, sum)

# Write raster
writeRaster(Total_cons, here::here("Data", "WaterUse", "Consumption", "TotalSum_2010.tif"), 
            format = "GTiff", overwrite = T)