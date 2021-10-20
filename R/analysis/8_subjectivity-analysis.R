################################################################################
# @Manuscript - "Hotspots of social and ecological impacts from freshwater stress and storage loss" (Huggins et al.) 
# @Description - Perform sensitivity analysis on the impact of subjective methodological decisions.
################################################################################

# Load 'here' library for easy path management.  
library(here)

# Import all setup and user-defined functions in R/setup and R/udfs folders
invisible(sapply(paste0(here("R/setup"), "/", list.files(here("R/setup"))), source)) 
invisible(sapply(paste0(here("R/udfs"), "/", list.files(here("R/udfs"))), source))

# Set methodology alternative labels
disL <- c('L3', 'L4', 'L5') # HydroBASINS levels 3, 4, or 5
sflw <- c('GSCD', 'GRUN') # GSCD or GRUN for streamflow
useC <- c('C', 'W') # Use consumption or withdrawal rates for stress calculations
agg <- c('Am', 'Fz') # Aggregate adaptive capacity and ecological sensitivity using fuzzy sum or arithmetic average

# Create tibble with all combinations of above alternatives
Combs <- expand.grid(disL, agg, useC, sflw) %>%  
  set_colnames(c('disL', 'agg', 'useC', 'sflw'))

# Flag zones of earthquake interference of TWS trends
basid <- raster(here('Data/HyBas4_cleaned.tif'))
eq_id <- c(4130, 4440, # Japan
           4451, 4430, 4457, 4457, 4459, 4460,
           5111, 5112, 5113, 5114, 5115, 5116, 5117, 5118, 5120, 5130) # Malay Pen.
eq_rm <- raster(WGS84_areaRaster(0.5)) # Initiate this raster

eq_rm[basid >= 0] <- 0

for (e in eq_id) {
  eq_rm[basid == e] <- 1
}

# Execute hotspot analysis for each alternative

for (z in 1:nrow(Combs)) {
  
  # Import data unaffected by alternatives ----
  tws  <- raster(here('Data/Rodell_GRACE_TWSt_raw.tif')) * 10 #cm to mm
  efn <- raster(here("Data/Dimensions/EFN_0d5.tif"))
  vsi <- raster(here("Data/Dimensions/VSI_aetw_0d5.tif"))
  ac15 <- raster(here::here("Data", "AdaptiveCap2015.tif"))
  
  
  # Import discretization scheme according to selection ----
  if (Combs$disL[z] == 'L3') {
    dis  <- raster(here('Data/HyBas3_cleaned.tif'))  
  }
  
  if (Combs$disL[z] == 'L4') {
    dis  <- raster(here('Data/HyBas4_cleaned.tif'))
  }
  
  if (Combs$disL[z] == 'L5') {
    dis  <- raster(here('Data/HyBas5_cleaned.tif')) 
  }
  
  # Import consumption or withdrawal rate according to selection ----
  if (Combs$useC[z] == 'W') {
    wuse <- raster(here('Data/Dimensions/TotalWithdrawls_2010.tif'))
  }
  
  if (Combs$useC[z] == 'C') {
    wuse <- raster(here('Data/WaterUse/Consumption/TotalSum_2010.tif')) 
  }
  
  # import streamflow dataset according to selection ----
  if (Combs$sflw[z] == 'GSCD') {
    roff <- raster(here('Data/GSCD_Qmean_0d5.tif'))
  }
  if (Combs$sflw[z] == 'GRUN') {
    roff <- raster(here('Data/GRUN/MeanGRUN.tif')) * 365.25 # /day to /year
  }
  
  # Generate inputs for the hotspot analysis, and ensure data from all alternatives can be executed without errors ----
  
  # Ensure data availability in all basins in whichever scheme being used
  ac.ext <- ac15
  ac.ext[ac.ext >= 0] <- 1
  ac.ext[is.na(ac15)] <- 0
  
  clip_df <- raster::stack(dis, ac.ext, WGS84_areaRaster(0.5)) %>% 
    as.data.frame() %>% 
    set_colnames(c('dis', 'ac.ext', 'area'))
  
  clip_df <- clip_df[complete.cases(clip_df$dis),]
  clip_df <- clip_df %>% 
    group_by(dis) %>% 
    summarize(
      cvg = weighted.mean(x = ac.ext, w = area, na.rm = T)
    )
  
  disclude <- clip_df %>% filter(cvg == 0) %>% pull(dis)
  
  if (length(disclude) > 0) {
    for (i in 1:length(disclude)) {
      dis[dis == disclude[i]] <- NA
    } 
  }
  
  # Establish a masking extent
  dis.ext <- dis
  dis.ext[dis.ext >= 0] <- 1
  
  # Derive ecological sensitivity inputs and indicator
  ef.ptl <- RasterAreaPercentiles(RasterToClassify = efn,
                                  WeightRaster = WGS84_areaRaster(0.5),
                                  MaskRaster = dis.ext,
                                  clipToExtent = "clip",
                                  CRS.set = crs(dis.ext),
                                  ext.set = extent(dis.ext))
  ef.ptl[ef.ptl > 0] <- (1 - ef.ptl[ef.ptl > 0]) + 0.01 # invert as high depths are less sensitive
  
  vsi[is.na(vsi) | dis.ext != 1] <- NA
  vsi.ptl <- RasterAreaPercentiles(RasterToClassify = vsi,
                                   WeightRaster = WGS84_areaRaster(0.5),
                                   MaskRaster = dis.ext,
                                   clipToExtent = "clip",
                                   CRS.set = crs(dis.ext),
                                   ext.set = extent(dis.ext))
  
  # Invert adaptive capacity
  ac <- 1 - ac15
  
  # Derive basin freshwater status inputs
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
  
  # Storage trend indicator
  tws_ind <- max(min((tws/(0.4*roff)), 1), -1)*-1
  
  # basin freshwater status:
  cind <- max(min(( (fws_ind + tws_ind)/2 ), 1), 0)
  
  
  # Set basin freshwater status to just the freshwater stress indicator where sufficient earthquake interference occurs
  eq_fix <- raster::stack(dis, eq_rm, WGS84_areaRaster(0.5)) %>% 
    as.data.frame() %>% 
    set_colnames(c('dis', 'eq_rm', 'area'))
  
  eq_fix <- eq_fix[complete.cases(clip_df$dis),]
  eq_fix <- eq_fix %>% 
    group_by(dis) %>% 
    summarize(
      eq_cov = weighted.mean(x = eq_rm, w = area, na.rm = T)
    )
  
  eq_ow <- eq_fix %>% filter(eq_cov > 0.1) %>% pull(dis)
  
  for (q in eq_ow) {
    cind[dis == q] <- fws_ind[dis == q]
  }
  
  # Crop to where the basin freshwater status exists
  cind.ext <- cind
  cind.ext[cind.ext >= 0] <- 1
  cind.ext[is.na(cind)] <- 0
  
  clip_df <- raster::stack(dis, cind.ext, WGS84_areaRaster(0.5)) %>% 
    as.data.frame() %>% 
    set_colnames(c('dis', 'cind.ext', 'area'))

  clip_df <- clip_df[complete.cases(clip_df$dis),]
  
  clip_df <- clip_df %>% 
    group_by(dis) %>% 
    summarize(
      cvg = weighted.mean(x = cind.ext, w = area, na.rm = T)
    )
  
  disclude <- clip_df %>% filter(cvg == 0) %>% pull(dis)
  
  if (length(disclude > 0)) {
    for (i in 1:length(disclude)) {
      dis[dis == disclude[i]] <- NA
    } 
  }
  
  # Convert raster stack to dataframe for faster computation
  c_df <- raster::stack(dis, cind, ef.ptl, vsi.ptl, ac, WGS84_areaRaster(0.5)) %>% 
    as.data.frame() %>% 
    set_colnames(c('dis', 'cind', 'efn', 'vsi', 'ac', 'area'))
  
  c_df <- c_df[complete.cases(c_df$dis),]
  
  c_df <- c_df %>% 
    group_by(dis) %>% 
    summarise(
      cind = weighted.quantile(x = cind, w = area, probs = 0.5, na.rm = T),
      efn = weighted.mean(x = efn, w = area, na.rm = T),
      vsi = weighted.mean(x = vsi, w = area, na.rm = T),
      ac = weighted.mean(x = ac, w = area, na.rm = T)
    )
  
  # Calculate ecological sensitivity
  c_df$ecosens <- (c_df$vsi + c_df$efn)/2
  c_df$ecosens <- c_df$ecosens/max(c_df$ecosens, na.rm = T)
  
  # choose here between fuzzy sum and arithmetic average for social-ecological sensitivity, according to selection ----
  if (Combs$agg[z] == 'Fz') {
    c_df$overallsens <- 1 - (1-c_df$ecosens)*(1-c_df$ac) # this is fuzzy sum
  }
  if (Combs$agg[z] == 'Am') {
    c_df$overallsens <- (c_df$ecosens+c_df$ac)/2 # this is mean
  }
  
  # Calculate social-ecological vulnerability
  c_df$overallprod <- c_df$cind * c_df$overallsens
  
  # Identify class breaks in vulnerability using Head/Tail scheme
  htb_o <- ht_breaks2.0(x = c_df$overallprod, tsh = 0.8)
  
  
  # Generate maps using these class breaks ----
  
  # Adaptive capacity and ecological sensitivity per basin
  ac.f <- FeatureAreaAverage(FT.id = dis, 
                             RawDS = ac,
                             AreaDS = WGS84_areaRaster(0.5), 
                             operation = 'mean',
                             varnam = 'ac')
  
  es.f <- FeatureAreaAverage(FT.id = dis, 
                             RawDS = (ef.ptl+vsi.ptl)/2,
                             AreaDS = WGS84_areaRaster(0.5), 
                             operation = 'mean',
                             varnam = 'ecosens')
  es.f <- es.f/max(es.f[], na.rm = T)
  
  # choose here between fuzzy sum and arithmetic average for social-ecological sensitivity, according to selection ----
  if (Combs$agg[z] == 'Fz') {
    ov.f <- 1 - (1 - ac.f)*(1-es.f) # fuzzy sum
  }
  if (Combs$agg[z] == 'Am') {
    ov.f <- (ac.f+es.f)/2 # mean 
  }
  
  # Classify each grid cell using Head/Tail breaks into transitional basins +, or not
  rclmat <- c(0, htb_o[1],        0, 
              htb_o[1], 1,        1) %>% 
    matrix(ncol = 3, byrow = T)
  ovhot <- reclassify(ov.f*cind, rclmat, include.lowest = T)
  
  # Write classified raster for future analysis
  writeRaster(ovhot, paste0(here('Data/Sensitivity'), "/", 
                            Combs$disL[z], "_", 
                            Combs$agg[z], "_", 
                            Combs$useC[z], "_", 
                            Combs$sflw[z], "_transitional", sep = ""), 
              format = 'GTiff', overwrite = T)
  
  # Classify each grid cell using Head/Tail breaks into hotspot basin or not
  rclmat <- c(0, htb_o[2],        0, 
              htb_o[2], 1,        1) %>% 
    matrix(ncol = 3, byrow = T)
  ovhot <- reclassify(ov.f*cind, rclmat, include.lowest = T)
  
  # Write classified raster for future analysis
  writeRaster(ovhot, paste0(here('Data/Sensitivity'), "/", 
                            Combs$disL[z], "_", 
                            Combs$agg[z], "_", 
                            Combs$useC[z], "_", 
                            Combs$sflw[z], "_hotspot", sep = ""),
              format = 'GTiff', overwrite = T)
  
  # Classify each grid cell into all four vulnerability classes, for web-app
  rclmat <- c(0,        htb_o[1],        0, 
              htb_o[1], htb_o[2],        1,
              htb_o[2], htb_o[3],        2,
              htb_o[3], Inf,             3) %>% 
    matrix(ncol = 3, byrow = T)
  VulClass <- reclassify(ov.f*cind, rclmat, include.lowest = T)
  
  VulClass_sf <- raster::rasterToPolygons(VulClass, na.rm = T, dissolve = T)
  
  # Write vulnerability class shapefile for web-app
  writeOGR(VulClass_sf, 
           dsn = here('Web-app'), 
           layer = paste0(Combs$disL[z], "_", 
                          Combs$agg[z], "_", 
                          Combs$useC[z], "_", 
                          Combs$sflw[z], sep = ""),
           driver = "ESRI Shapefile", overwrite = T)
  
  
  print(z)
  
} # End of loop


# Generate plots of the frequency of transitional and hotspot identification across the 24 methodology alternatives ---- 

# List all rater files that identify grid cells as transitional basins and hotspots
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

  tmap_save(tm, paste0(here('Data/Sensitivity/'), "/Freq", i,".svg", sep = ""),
            units = "in")
}