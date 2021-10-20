################################################################################
# @Manuscript - "Hotspots of social and ecological impacts from freshwater stress and storage loss" (Huggins et al.) 
# @Description - Derive the basin freshwater status indicator. 
################################################################################

# Load 'here' library for easy path management.  
library(here)

# Import all setup and user-defined functions in R/setup and R/udfs folders
invisible(sapply(paste0(here("R/setup"), "/", list.files(here("R/setup"))), source)) 
invisible(sapply(paste0(here("R/udfs"), "/", list.files(here("R/udfs"))), source))

# Import data ----
wuse <- raster(here('Data/Dimensions/TotalWithdrawls_2010.tif'))
roff <- raster(here('Data/GSCD_Qmean_0d5.tif'))
tws  <- raster(here('Data/Rodell_GRACE_TWSt_raw.tif')) * 10 # cm to mm
dis  <- raster(here('Data/HyBas4_cleaned.tif'))

# Calculate area-weighted basin average of all datasets
wuse <- FeatureAreaAverage(FT.id = dis, 
                           RawDS = wuse,
                           AreaDS = WGS84_areaRaster(0.5), 
                           operation = 'mean',
                           varnam = 'wuse')

roff <- FeatureAreaAverage(FT.id = dis, 
                           RawDS = roff,
                           AreaDS = WGS84_areaRaster(0.5), 
                           operation = 'mean', 
                           varnam = 'roff')

tws <- FeatureAreaAverage(FT.id = dis, 
                          RawDS = tws, 
                          AreaDS = WGS84_areaRaster(0.5), 
                          operation = 'mean', 
                          varnam = 'tws')

# Freshwater stress indicator
fws_ind <- min((wuse/(0.4*roff)), 1)

# Terrestrial water storage trend indicator
tws_ind <- max(min((tws/(0.4*roff)), 1), -1) * -1

# Basin freshwater status
cmb_ind <- max(min(( (fws_ind + tws_ind)/2 ), 1), 0)

# Earthquake interference affects TWS trends in Japan & Malay Peninsula. We set the combined indicator score to just the FW stress indicator score for these regions.
eqbas <- c(4130, 4440, # Japan
           4451, 4430, 4457, 4457, 4459, 4460,
           5111, 5112, 5113, 5114, 5115, 5116, 5117, 5118, 5120, 5130) # Malay Pen.

for (i in eqbas) {
  cmb_ind[dis == i] <- fws_ind[dis == i]
}

# Write the basin freshwater status indicator raster for future use           
writeRaster(cmb_ind, here('Data/fwss_ind_comb.tif'), format = 'GTiff',
            overwrite = T)

# Map the two inputs and the combined indicator ----

# Import vector boundary files
disborders <- sf::read_sf(here('Data/hybas_l4_0d5.shp'))
coastlines <- sf::read_sf(here('Data/hybas_l4_coastlines.shp'))

# Stack the rasters to loop through plotting function
plotstack <- stack(fws_ind, tws_ind, cmb_ind)

for (i in 1:nlayers(plotstack)) {
  
  # Reproject raster
  plt.obj <- tmap_clipproj(plotstack[[i]])
  
  # Specify color palette per indicator, upper and lower limits, and midpoint
  if (i == 1) {
    cus.pal <- scico(20, palette = "bilbao", direction = 1)
    ul <- 1; ll <- 0; mp <- 0.5
  } 
  if (i == 2) {
    cus.pal <- scico(20, palette = "vikO", direction = 1)[4:17]
    ul <- 1; ll <- -1; mp <- 0
    
  } 
  
  tm <- tm_shape(plt.obj, projection = "+proj=robin") +
    tm_raster(style = "cont", palette = cus.pal, 
              midpoint = mp, breaks = c(ll, ul), alpha = 1.00) +
    tm_shape(disborders) +
    tm_borders(lwd = 0.1, col = "grey40") +
    tm_shape(coastlines) +
    tm_borders(lwd = 0.7, col = "black") +
    tm_layout(legend.show = F, earth.boundary = c(-179, -60, 179, 88),
              earth.boundary.color = "white", space.color = "white",
              legend.frame = F, frame = F,
              outer.margins = c(-0, -0.09, -0, -0.04)) # B, L, T, R
  
  # Add a grey mask for low combined indicator values
  if (i == 3) {
    cus.pal <- scico(20, palette = "lajolla", direction = 1)[4:17]
    ul <- 1; ll <- 0.0; mp <- 0.5
    
    greymask <- plt.obj
    greymask[plt.obj >= 0.05] <- NA
    greymask[greymask >= 0] <- 1
    
    tm <- tm_shape(plt.obj, projection = "+proj=robin") +
      tm_raster(style = "cont", palette = cus.pal, 
                midpoint = mp, breaks = c(ll, ul), alpha = 1.00) +
      tm_shape(greymask) +    
      tm_raster(style = "cat", palette = "grey80") +
      tm_shape(disborders) +
      tm_borders(lwd = 0.1, col = "grey40") +
      tm_shape(coastlines) +
      tm_borders(lwd = 0.7, col = "black") +
      tm_layout(legend.show = F, earth.boundary = c(-179, -60, 179, 88),
                earth.boundary.color = "white", space.color = "white",
                legend.frame = F, frame = F,
                outer.margins = c(-0, -0.09, -0, -0.04)) # B, L, T, R
  }
  
  tm
  
  assign(paste0("plt_", i, sep = ""), tm)
  
  tmap_save(tm, 
            paste0("C:/Users/xande/Desktop/jt_prep/2_ind_step", i, ".svg", sep = ""),
            units = "in")
}