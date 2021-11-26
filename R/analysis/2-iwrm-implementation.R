# Name: 2-iwrm-implementation.R
# Description: Evaluate IWRM implementation levels per basin to enable comparison to vulnerability results.

# Use 'here' package for easy path management
library(here)

# Import all setup and user-defined functions in R/setup and R/udfs folders
invisible(sapply(paste0(here("R/setup"), "/", list.files(here("R/setup"))), source)) 
invisible(sapply(paste0(here("R/udf"), "/", list.files(here("R/udf"))), source))

# Import data ----
# Geospatial areal units of analysis
dis  <- raster(here('Data/HyBas4_cleaned.tif')) # basin scheme

# IWRM data
SDG651_17 <- readr::read_csv(here('Data/SDG651_2017_baseline_cleaned.csv'))
SDG651_20 <- readr::read_csv(here('Data/SDG651_2020_stripped.csv'))

# Country boundaries to merge IWRM data to
faoGAUL <- read_sf(here('Data/go_countryborders_faogaul2014Polygon.shp'))

# Joining IWRM data with w/ GAUL and write to shapefile ----
SDGgaul_17 <- merge(x = faoGAUL, y = SDG651_17, by.x = 'ccode', by.y = 'ISOcode')
SDGgaul_20 <- merge(x = faoGAUL, y = SDG651_20, by.x = 'ccode', by.y = 'ISO3')

write_sf(obj = SDGgaul_17, dsn = here('Data/'), layer = 'SDGgaul_17.shp', 
         driver = 'ESRI Shapefile')

write_sf(obj = SDGgaul_20, dsn = here('Data/'), layer = 'SDGgaul_20.shp', 
         driver = 'ESRI Shapefile')

# Rasterize GAUL dataset for transboundary calculations  ----
faoGAUL$uid <- seq(1, nrow(faoGAUL), 1)
write_sf(obj = SDGgaul_20, dsn = here('Data/'), layer = 'gaul_ids.shp', 
         driver = 'ESRI Shapefile')

gdalUtils::gdal_rasterize(src_datasource = here::here('Data/gaul_ids.shp'),
                          dst_filename = here::here('Data/gaul_ids.tif'),
                          a = 'uid',
                          te = c(-180, -90, 180, 90),
                          tr = c(0.05, 0.05),
                          a_nodata = NA)

gdalUtils::gdalwarp(srcfile = here::here('Data/gaul_ids.tif'),
                    dstfile = here::here('Data/gaul_ids_halfdegree.tif'),
                    te = c(-180, -90, 180, 90),
                    tr = c(0.5, 0.5),
                    r = 'mode')

# Rasterize IWRM implementation levels for both years of reporting ----
gdalUtils::gdal_rasterize(src_datasource = here::here('Data/SDGgaul_17.shp'),
                          dst_filename = here::here('Data/SDGgaul_17.tif'),
                          a = 'Final',
                          te = c(-180, -90, 180, 90),
                          tr = c(0.05, 0.05),
                          a_nodata = NA)

gdalUtils::gdal_rasterize(src_datasource = here::here('Data/SDGgaul_20.shp'),
                          dst_filename = here::here('Data/SDGgaul_20.tif'),
                          a = 'SDG651_',
                          te = c(-180, -90, 180, 90),
                          tr = c(0.05, 0.05),
                          a_nodata = NA)

# Resample to half degree using mean value
gdalUtils::gdalwarp(srcfile = here::here('Data/SDGgaul_17.tif'),
                    dstfile = here::here('Data/SDGgaul_17_halfdegree.tif'),
                    te = c(-180, -90, 180, 90),
                    tr = c(0.5, 0.5),
                    r = 'average')

gdalUtils::gdalwarp(srcfile = here::here('Data/SDGgaul_20.tif'),
                    dstfile = here::here('Data/SDGgaul_20_halfdegree.tif'),
                    te = c(-180, -90, 180, 90),
                    tr = c(0.5, 0.5),
                    r = 'average')

# Replace NAs in 2020 with values in 2017 if they exist
iwrm_2020 <- raster(here::here('Data/SDGgaul_20_halfdegree.tif'))
iwrm_2017 <- raster(here::here('Data/SDGgaul_17_halfdegree.tif'))
iwrm_fill <- iwrm_2020
iwrm_fill[is.na(iwrm_2020)] <- iwrm_2017[is.na(iwrm_2020)]
writeRaster(iwrm_fill, here::here('Data/iwrm_fill.tif'), overwrite = T)

# Determine IWRM implementation per basin
c_df <- raster::stack(dis, # basins
                      iwrm_fill, # IWRM data
                      raster(here::here('Data/gaul_id_halfdegree.tif')), # country ids 
                      WGS84_areaRaster(0.5)) %>% 
  as.data.frame() %>% 
  set_colnames(c('id', 'iwrm', 'natid','area'))
c_df <- c_df[complete.cases(c_df$id),]

c_df <- c_df %>% 
  group_by(id) %>%
  summarize(
    iwrm  = weighted.mean(x = iwrm, w = area, na.rm = T),
    tsbd = length(unique(natid))
)

# Import basins_data.shp and add iwrm to the attribute table
basins <- sf::read_sf(here('R/plotting+stats/Basin_data.shp'))

basins <- merge(x = basins, y = c_df, by.x = 'pfaf_id', by.y = 'id', all.x = T)

write_sf(obj = basins, dsn = here('R/plotting+stats/'), layer = 'Basin_iwrm_data.shp', 
         driver = 'ESRI Shapefile')