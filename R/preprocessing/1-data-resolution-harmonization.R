# Name: 1-data-resolution-harmonization.R
# Description: Harmonize input data to 0.5° resolution 

# Use 'here' package for easy path management
library(here)

# Import all setup scripts and user-defined functions in R/setup and R/udf folders
invisible(sapply(paste0(here("R/setup"), "/", list.files(here("R/setup"))), source)) 
invisible(sapply(paste0(here("R/udf"), "/", list.files(here("R/udf"))), source))

# Notes:
# Freshwater consumption & withdrawal rates and population data are directly available in 0.5° resolution
# IWRM data is handled in script 5-iwrm-implementation-comparison.R
# Streamflow data are prepared in Preprocessing/3-streamflow-alternatives.R

# TWS trends ----
Rodell_raw <- readr::read_csv(here::here("Data", "41586_2018_123_MOESM1_ESM.csv"),
                              col_names = F) %>%
  as.matrix()
TWSt.ras <- raster(Rodell_raw)
extent(TWSt.ras) <- c(-180, 180, -90, 90)
crs(TWSt.ras) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"
TWSt.ras <- flip(TWSt.ras, direction = "y")
writeRaster(TWSt.ras, here::here("Data", "Rodell_GRACE_TWSt_raw.tif"),
            format = "GTiff", overwrite = T)

# Vegetation sensitivity index - water availability ----
gdalwarp(srcfile = "D:/!! Geodatabase/VSI/SensAETW.tif",
         dstfile = here::here("Data", "Dimensions", "VSI_aetw_0d5.tif"),
         te = c(-180, -90, 180, 90),
         tr = c(0.5, 0.5),
         r = "average",
         output_Raster = TRUE,
         overwrite = TRUE,
         verbose = TRUE)                                     

# Environmental flow limit sensitivity to groundwater head decline ----
EFN.dg <- raster(here::here("Data", "Dimensions", "headdrop2limit_hydr06.map"))
x = 0.5/res(EFN.dg)
EFN.dg_0d5 <- raster::aggregate(x = EFN.dg, fact = x, fun = mean, expand = F, na.rm = T, 
                                  filename = here::here("Data", "Dimensions", 
                                                        "EFN_0d5.tif"),
                                overwrite = T)

# Adaptive capacity ----
AdaptiveCap <- nc_open(here::here("Data/adaptive_capacity.nc"))
lon <- ncvar_get(AdaptiveCap, "longitude")
lat <- ncvar_get(AdaptiveCap, "latitude")
dname <- "adaptive_capacity"
ncatt_get(AdaptiveCap, dname, "long_name")
ncatt_get(AdaptiveCap, dname, "units")
AdaptiveCap <- ncvar_get(AdaptiveCap, dname)
dim(AdaptiveCap)
AdaptiveCap_2015 <- raster(t(AdaptiveCap[,,26]), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), 
                           crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
extent(AdaptiveCap_2015) <- c(-180, 180, -90, 90)

x = 0.5/res(AdaptiveCap_2015)
AC_2015_0d5 <- raster::aggregate(x = AdaptiveCap_2015, fact = x, fun = mean, expand = F, na.rm = T, 
                                filename = here::here("Data", "AdaptiveCap2015.tif"),
                                overwrite = T)

# Food crop calories ----
Cal <- raster("D:/!! Geodatabase/Agriculture/SPAM/SPAM_2010_foodcrops_kcal.tif")
x = 0.5/res(Cal)
Cal_0d5 <- raster::aggregate(x = Cal, fact = x, fun = sum, expand = F, na.rm = T, 
                             filename = here::here("Data", "Dimensions", "SPAM_foodcrops_kcal_0d5.tif"))

# GDP ----
GDP <- nc_open(here::here("Data", "Dimensions", "GDP_PPP_1990_2015_5arcmin_v2.nc"))
lon = ncvar_get(GDP, "longitude")
lat = ncvar_get(GDP, "latitude")
v.name = "GDP_PPP"
ncatt_get(GDP, v.name, "long_name")
ncatt_get(GDP, v.name, "units")
fv = ncatt_get(GDP, v.name, "_FillValue")
GDP_years <- ncvar_get(GDP, v.name)
dim(GDP_years)

# 2015 is year 26 (last year) in netCDF
GDP_2015 <- raster(t(GDP_years[,,26]), xmn=min(lon), xmx=max(lon), 
                   ymn=min(lat), ymx=max(lat), 
                   crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
extent(GDP_2015) <- c(-180, 180, -90, 90)
x = 0.5/res(GDP_2015)
GDP_2015_0d5 <- raster::aggregate(x = GDP_2015, fact = x, fun = sum, expand = F, na.rm = T, 
                                  filename = here::here("Data", "Dimensions", 
                                                        "GDP_2015_0d5.tif"))


# Amphibian species richness ----
amph <- raster(here::here("Data", "Dimensions", "all_amphibians.tif"))
extent(amph) = c(-180, 180, -90, 90)
x = 0.5/res(amph)
# Use maximum taxonomic biodiversity value per larger grid cell. Using gdalwarp instead of raster::resample function as more efficient.
gdalwarp(srcfile = here::here("Data", "Dimensions", "all_amphibians.tif"),
         dstfile = here::here("Data", "Dimensions", "amph_0d5.tif"),
         te = c(-180, -90, 180, 90),
         tr = c(0.5, 0.5),
         r = "max",
         output_Raster = TRUE,
         overwrite = TRUE,
         verbose = TRUE)


# Ramsar sites ----
ramsar <- readr::read_csv('D:/!! Geodatabase/Biodiversity/ris-tabular-20210322.csv')
ramsar.pts <- SpatialPoints(cbind(ramsar$Longitude, ramsar$Latitude))
crs(ramsar.pts) <- crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+")

ramsar.pts <- raster::rasterize(ramsar.pts, WGS84_areaRaster(0.5), 
                               fun = 'count')

writeRaster(ramsar.pts, here("Data/ramsar_count.tif"),
            format = 'GTiff', overwrite = T)


# World Conflict Chronology ---- 
WCC <- readr::read_csv('D:/!! Geodatabase/Social-data/Water-conflict/WCC/WCC_PacificInstitute_926conflicts.csv') %>% 
  as.data.frame() %>%
  filter(Start >= 2000)
WCC$Long %<>% as.numeric()
WCC <- WCC[complete.cases(WCC$Lat),]
WCC <- WCC[complete.cases(WCC$Long),]
WCC <- st_as_sf(WCC, coords = c("Long", "Lat"), crs = 4326)
WCC$Counter <- 1

WCC_L20 <- raster::rasterize(WCC, WGS84_areaRaster(0.5), 
                                   field = 'Counter', fun = 'count')
writeRaster(WCC_L20, 
            here::here("Data", "Dimensions", "WCC_since2000.tif"),
            format = 'GTiff', overwrite = T)