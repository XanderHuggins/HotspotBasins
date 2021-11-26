# Name: 4-hydrobasins-grid-snap-and-clean.R
# Description: Preprocess HydroBASINS for use by reconciling to 0.5° grid and cleaning for sufficient data coverage. At the bottom of this script, global statistics are calculated. 

# Use 'here' package for easy path management
library(here)

# Import all setup and user-defined functions in R/setup and R/udfs folders
invisible(sapply(paste0(here("R/setup"), "/", list.files(here("R/setup"))), source)) 
invisible(sapply(paste0(here("R/udf"), "/", list.files(here("R/udf"))), source))

# NOTE: Only the preprocessing of Level 4 is shown below. However, it is repeated for Levels 3 and 5 for use in the sensitivity analysis in 'analysis/5-sensitivity-to-methods.R'

# Import all world regions
af <- sf::read_sf('D:/!! Geodatabase/Basins/hybas/Raw/hybas_af_lev04_v1c.shp')
ar <- sf::read_sf('D:/!! Geodatabase/Basins/hybas/Raw/hybas_ar_lev04_v1c.shp')
as <- sf::read_sf('D:/!! Geodatabase/Basins/hybas/Raw/hybas_as_lev04_v1c.shp')
au <- sf::read_sf('D:/!! Geodatabase/Basins/hybas/Raw/hybas_au_lev04_v1c.shp')
eu <- sf::read_sf('D:/!! Geodatabase/Basins/hybas/Raw/hybas_eu_lev04_v1c.shp')
gr <- sf::read_sf('D:/!! Geodatabase/Basins/hybas/Raw/hybas_gr_lev04_v1c.shp')
na <- sf::read_sf('D:/!! Geodatabase/Basins/hybas/Raw/hybas_na_lev04_v1c.shp')
sa <- sf::read_sf('D:/!! Geodatabase/Basins/hybas/Raw/hybas_sa_lev04_v1c.shp')
si <- sf::read_sf('D:/!! Geodatabase/Basins/hybas/Raw/hybas_si_lev04_v1c.shp')

# Merge into single file
hybas_list <- list(af, ar, as, au, eu, gr, na, sa, si)
hybas_l_set   <- do.call(rbind, hybas_list) 

# Write combined shapefile
write_sf(obj = hybas_l_set, dsn = here('Data/'), layer = 'hybas_l4.shp', 
         driver = 'ESRI Shapefile')

# Rasterize, using finer resolution then aggregating to avoid clipping
gdalUtils::gdal_rasterize(src_datasource = here::here('Data/hybas_l4.shp'),
                          dst_filename = here::here('Data/hybas_l4_raw.tif'),
                          a = 'PFAF_ID',
                          te = c(-180, -90, 180, 90),
                          tr = c(0.05, 0.05),
                          a_nodata = NA)

# Aggregate to half degree using modal value
gdalUtils::gdalwarp(srcfile = here::here('Data/hybas_l4_raw.tif'),
                    dstfile = here::here('Data/hybas_l4_pfaf_halfdegree.tif'),
                    te = c(-180, -90, 180, 90),
                    tr = c(0.5, 0.5),
                    r = 'mode')
                    
# Convert raster to shapefile, for plotting
hybas_0d5 <- rasterToPolygons(raster(here::here('Data/hybas_l4_pfaf_halfdegree.tif')), 
                              dissolve = T)
hybas_0d5 <- sf::st_as_sf(hybas_0d5)
write_sf(obj = hybas_0d5, dsn = here('Data/'), layer = 'hybas_l4_0d5.shp', 
         driver = 'ESRI Shapefile')

# Clean basins for analysis, by masking regions with partial input data coverage, removing Greenland
hybas_rasag <- raster(here::here('Data/hybas_l4_pfaf_halfdegree.tif'))

# GRUN has limited data coverage in high latitudes, so use as extent limiting dataset
GRUN <- raster(here("Data", "GRUN", "MeanGRUN.tif"))

# Create binary raster indicating where GRUN exists
GRUNexist <- raster(GRUN)
GRUNexist[] <- 0
GRUNexist[GRUN >= 0] <- 1

# Calculate GRUN coverage per basin
GRUN_coverage <- FeatureAreaAverage(FT.id = hybas_rasag,
                                    RawDS = GRUNexist,
                                    AreaDS = WGS84_areaRaster(0.5),
                                    operation = "mean",
                                    varnam = "GRUNcoverage")

# Exclude basins with low GRUN coverage
GRUN_clip <- GRUN_coverage
GRUN_clip[GRUN_clip < 0.67] <- NA
GRUN_clip[GRUN_clip > 0] <- 1
plot(GRUN_clip)

# Manually remove other problematic basins

# All of Greenland
manualclip <- raster(hybas_rasag)
manualclip[] <- 1
manualclip[hybas_rasag >= 9000] <- NA # Greenland PFAF ids (need to change this to 900 for level 3, and 90000 for level 5)

# Combine clipping extents
ClipExtent <- GRUN_clip * manualclip

# Clip the HydroBASIN ID raster
hybas_rasag_cleaned <- hybas_rasag*ClipExtent

# Basins to exclude because of insufficient adaptive capacity coverage (manually identified)
# Only done for Level 4
exclude <- readr::read_csv(file = here::here("Data/basin_lists_1.csv"),
                            col_select = 'insufficient_AC_coverage',
                            col_names = T)

for (i in exclude$insufficient_AC_coverage) {
  hybas_rasag_cleaned[hybas_rasag_cleaned == i] <- NA
}

# Calculate how many basins are retained in this final version
hybas_rasag_cleaned %>% unique() %>% length()
## output = 1204

# Write cleaned raster
writeRaster(hybas_rasag_cleaned, here::here("Data", "HyBas4_cleaned.tif"), 
            format = "GTiff", overwrite = T)


# Calculate coverage statistics for Level 4 HydroBASINS
hbl4 <- sf::read_sf("D:/!! Geodatabase/Basins/hybas/Processed/hybas_l4.shp")

# General stats
mean(hbl4$SUB_AREA)
median(hbl4$SUB_AREA)

# How many basins in original set
hbl4 %>% pull(PFAF_ID) %>% unique() %>% length() 
## output = 1341 (number of basins in original set)

# What is surface area of full set
c_df <- raster::stack(here::here('Data/hybas_l4_pfaf_halfdegree.tif'),
                      WGS84_areaRaster(0.5)) %>% 
  as.data.frame() %>% 
  set_colnames(c('id', 'area'))

c_df %>% filter(id > 0) %>% pull(area) %>% sum()/1e6
## output = 148.2 (surface area of original set estimated after rasterization)

# Now compare to cleaned basin template used in study
l4_cleaned_ras <- raster(here("Data/HyBas4_cleaned.tif"))

c_df <- raster::stack(l4_cleaned_ras, WGS84_areaRaster(0.5)) %>% 
  as.data.frame() %>% 
  set_colnames(c('id', 'area'))

temp <- c_df %>% 
  group_by(id) %>% 
  summarize(
    sa = sum(area, na.rm = T)
    ) 

temp <- temp[complete.cases(temp$id),]
summary(temp$sa) # Calculates average basin size, median = 71,391 km2
nrow(temp) # Calculates number of basins, output = 1204 

# Basin retention:
1204/1341
## output = 89.8%

# Calculate area retention
c_df %>% filter(id > 0) %>% pull(area) %>% sum()/1e6
## output = 143.2
143.2/148.2
# output = 96.6%
