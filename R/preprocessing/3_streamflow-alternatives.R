################################################################################
# @Manuscript - "Hotspots of social and ecological impacts from freshwater stress and storage loss" (Huggins et al.) 
# @Description - Generate annual streamflow raster from (1) GSCD and (2) GRUN datasets.
################################################################################

# Load 'here' library for easy path management.  
library(here)

# Import all setup and user-defined functions in R/setup and R/udfs folders
invisible(sapply(paste0(here("R/setup"), "/", list.files(here("R/setup"))), source)) 
invisible(sapply(paste0(here("R/udfs"), "/", list.files(here("R/udfs"))), source))

# Convert GSCD netcdf to raster ----
GSDCncdf <- nc_open('D://!! Geodatabase/GSDC/GSCD_v2.0.nc')
lon <- ncvar_get(GSDCncdf, "longitude")
lat <- ncvar_get(GSDCncdf, "latitude")
dname <- "QMEAN"
gsdc_QMEAN <- ncvar_get(GSDCncdf, dname)
gsdc_QMEAN <- raster(gsdc_QMEAN, xmn=min(lon), xmx=max(lon), 
                     ymn=min(lat), ymx=max(lat), 
                     crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ 
                           towgs84=0,0,0"))
extent(gsdc_QMEAN) <- extent(c(-180, 180, -90, 90))

# As GSCD is provided at 0.05 degrees, to we aggregate to 0.5 degrees but weight contributing cell by cell size for more accurate aggregations
cellID <- WGS84_areaRaster(0.5)
cellID[] <- seq(1, ncell(cellID[]), by = 1)
cellID0d05 <- raster::resample(x = cellID, y = WGS84_areaRaster(0.05), method = 'ngb')

c.df <- data.frame(cellID0d05[], gsdc_QMEAN[], WGS84_areaRaster(0.05)[]) %>% 
  set_colnames(c('id', 'Qmean', 'area'))
c.df <- c.df[complete.cases(c.df$dis),]

c.df <- c.df %>% 
  group_by(id) %>% 
  summarize(
    Qmean = weighted.mean(x = Qmean, w = area, na.rm = T)
  )

# format for raster reclassification
c.df$uplim <- c.df$id+0.5
c.df$lowlim <- c.df$id-0.5
rclmtx <- matrix(c(c.df$lowlim, c.df$uplim, c.df$Qmean), ncol = 3)

Qmean_0d5 <- reclassify(cellID, rclmtx)
plot(Qmean_0d5)

writeRaster(Qmean_0d5, here('Data/GSCD_Qmean_0d5.tif'),
            format = 'GTiff', overwrite = T)


# Calculate GRUN annual streamflow over the period 2000-2010 ----

# Initiate area raster
WGSarea <- WGS84_areaRaster(0.5)

# Import GRUN netcdf file 
GRUN <- nc_open("D:/!! Geodatabase/GRUN/GRUN_v1_GSWP3_WGS84_05_1902_2014.nc")
lon <- ncvar_get(GRUN, 'X')
lat <- ncvar_get(GRUN, 'Y')
vname <- 'Runoff'
fillval <- ncatt_get(GRUN, vname, '_FillValue')
GRUN_ext <- ncvar_get(GRUN, vname)
GRUN_ext[GRUN_ext == fillval$value] <- NA # Set fill values to NA
dim(GRUN_ext)[3]/(12*(2015-1902)) # Check right no. months, from start 1902 to end 2014

# Write raster for sum of gridded runoff for each year 2000-2010
M0_2000 <- 12*(1999-1902) + 1 # identify index of first month of 2000
DperM <- c(31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

for (i in 1:10) {
  
  j = (M0_2000) + (12*(1-i))
  AnnSum <- raster()
  # extract & sum each month in calendar year
  for (k in 1:12) {
    if (k == 1) {
      AnnSum <- DperM[k] * raster(t(GRUN_ext[,,(j+k-1)]), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
    }
    AnnSum <- AnnSum + 
      (DperM[k] * raster(t(GRUN_ext[,,(j+k-1)]), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0")))
  }
  
  # average out to mm/d 
  AnnSum <- AnnSum/sum(DperM)
  
  # flip y axis
  AnnSum <- flip(AnnSum, direction = 'y')
  
  # fill extent
  extent(AnnSum) <- c(-180, 180, -90, 90)
  
  # write each year as .tif file
  writeRaster(AnnSum, paste(here::here("Data", "GRUN"), "/", 1999+i, "_totalrunoff.tif", sep=""), 
              format = "GTiff",overwrite = T)
  
}

# Import all written rasters, and calculate mean value per grid over 10-year window
temp <- list.files(path = here("Data", "GRUN"), pattern = '_totalrunoff.tif$', full.names = T)
GRUNyears <- lapply(temp, raster) %>% raster::stack()
MeanGRUN <- calc(GRUNyears, fun = mean)
names(MeanGRUN) <- "meanGRUN_2000_2010"

# Write raster
writeRaster(MeanGRUN, here::here("Data", "GRUN", "MeanGRUN.tif"), format = "GTiff", overwrite = T)
