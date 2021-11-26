# Name: figure-S10-plotting.R
# Description: Generate plots of the frequency of identification as transitional and hotspot basins across 24 methodology alternatives 

# Use 'here' package for easy path management
library(here)

# Import all setup and user-defined functions in R/setup and R/udfs folders
invisible(sapply(paste0(here("R/setup"), "/", list.files(here("R/setup"))), source)) 
invisible(sapply(paste0(here("R/udf"), "/", list.files(here("R/udf"))), source))

# List all raster files that identify grid cells as transitional basins and hotspots
imp <- list.files(path = here('Data/Sensitivity/'), pattern = '_transitional')

# Initialize frequency raster
FreqTrans <- raster(WGS84_areaRaster(0.5))
FreqTrans[] <- 0

# Sum the frequencies across all rasters
for (i in 1:length(imp)) {
  FreqTrans <- FreqTrans + raster(paste0(here('Data/Sensitivity/'), 
                                         "/", imp[i], sep = ""))
}

# List all rater files that identify grid cells as transitional basins and hotspots
imp <- list.files(path = here('Data/Sensitivity/'), pattern = '_hotspot')

# Initialize frequency raster
FreqHot <- raster(WGS84_areaRaster(0.5))
FreqHot[] <- 0

# Sum the frequencies across all rasters
for (i in 1:length(imp)) {
  FreqHot <- FreqHot + raster(paste0(here('Data/Sensitivity/'), 
                                     "/", imp[i], sep = ""))
}

# Plot maps ----

# Import vector boundary files
coastlines <- sf::read_sf(here('Data/hybas_l4_coastlines.shp'))

# Loop through both maps
plotstack <- stack(FreqTrans, FreqHot)

for (i in 1:nlayers(plotstack)) {
  # Reproject for plotting in Robinson projection
  plt.obj <- tmap_clipproj(plotstack[[i]])
  
  # Color palette
  cus.pal <- c('grey80', scico(n = 25, alpha = 0.8, 
                               palette = 'batlow', dir = -1)[2:25])
  
  # Emphasize basins with outline
  emph <- raster(plotstack[[i]])
  for (j in c(seq(1, 24, by = 1))) 
  {emph[plotstack[[i]] == j] <- 1}
  
  emph <- raster::rasterToPolygons(emph, na.rm = T, dissolve = T)
  
  # Plot map
  tm <- tm_shape(plt.obj, projection = "+proj=robin") + tm_raster(style = "cat", palette = cus.pal) +
    tm_shape(coastlines) + tm_borders(lwd = 0.7, col = "black") +
    tm_shape(emph) + tm_borders(lwd = 1, col = "black") +
    tm_layout(legend.show = F, earth.boundary = c(-179, -60, 179, 88),
              earth.boundary.color = "white", space.color = "white",
              legend.frame = F, frame = F,
              outer.margins = c(-0, -0.09, -0, -0.04)) # B, L, T, R
  tm
  
  tmap_save(tm, paste0(here('plot_outputs/'), "/Freq", i,".svg", sep = ""),
            units = "in")
}
