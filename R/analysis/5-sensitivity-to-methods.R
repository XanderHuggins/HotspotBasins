# Name: 5-sensitivity-to-methods.R
# Description: Reproduce hotspot basin results for 24 methodological alternatives to consider the sensitivity of hotspot basin results to methodological decisions

# Use 'here' package for easy path management
library(here)

# Import all setup and user-defined functions in R/setup and R/udfs folders
invisible(sapply(paste0(here("R/setup"), "/", list.files(here("R/setup"))), source)) 
invisible(sapply(paste0(here("R/udf"), "/", list.files(here("R/udf"))), source))

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


eq_intfr <- readr::read_csv(file = here::here("Data/basin_lists_1.csv"),
                            col_select = 'eq_interference',
                            col_names = T,
                            show_col_types = F)
eq_rm <- raster(WGS84_areaRaster(0.5)) # Initiate this raster

eq_rm[basid >= 0] <- 0

for (r in eq_intfr$eq_interference) {
  eq_rm[basid == r] <- 1
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
  fw_status <- max(min(( (fws_ind + tws_ind)/2 ), 1), 0)
  
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
    fw_status[dis == q] <- fws_ind[dis == q]
  }
  
  # Crop to where the basin freshwater status exists
  fw_status.ext <- fw_status
  fw_status.ext[fw_status.ext >= 0] <- 1
  fw_status.ext[is.na(fw_status)] <- 0
  
  clip_df <- raster::stack(dis, fw_status.ext, WGS84_areaRaster(0.5)) %>% 
    as.data.frame() %>% 
    set_colnames(c('dis', 'fw_status.ext', 'area'))
  
  clip_df <- clip_df[complete.cases(clip_df$dis),]
  
  clip_df <- clip_df %>% 
    group_by(dis) %>% 
    summarize(
      cvg = sum(fw_status.ext, na.rm = T)
    )
  
  disclude <- clip_df %>% filter(cvg == 0) %>% pull(dis)
  
  if (length(disclude > 0)) {
    for (i in 1:length(disclude)) {
      dis[dis == disclude[i]] <- NA
    } 
  }
  
  # Convert raster stack to dataframe for faster computation
  c_df <- raster::stack(dis, fw_status, ef.ptl, vsi.ptl, ac15, WGS84_areaRaster(0.5)) %>% 
    as.data.frame() %>% 
    set_colnames(c('dis', 'fw_status', 'efn', 'vsi', 'ac', 'area'))
  
  c_df <- c_df[complete.cases(c_df$dis),]
  
  c_df <- c_df %>% 
    group_by(dis) %>% 
    summarise(
      fw_status = weighted.quantile(x = fw_status, w = area, probs = 0.5, na.rm = T),
      efn = weighted.mean(x = efn, w = area, na.rm = T),
      vsi = weighted.mean(x = vsi, w = area, na.rm = T),
      ac = weighted.mean(x = ac, w = area, na.rm = T)
    )
  
  # Calculate ecological sensitivity
  c_df$ecol_sens <- (c_df$vsi + c_df$efn)/2
  c_df$ecol_sens <- c_df$ecol_sens/max(c_df$ecol_sens, na.rm = T)
  
  # Show that inverted adaptive capacity is used as social sensitivity indicator
  c_df$socl_sens <- 1 - c_df$ac
  
  # choose here between fuzzy sum and arithmetic average for social-ecological sensitivity, according to selection ----
  if (Combs$agg[z] == 'Fz') {
    c_df$ses_sens <- 1 - ((1 - c_df$ecol_sens) * (1- c_df$socl_sens)) # this is fuzzy sum
  }
  if (Combs$agg[z] == 'Am') {
    c_df$ses_sens <- (c_df$ecol_sens + c_df$socl_sens)/2 # this is mean
  }
  
  # Calculate social-ecological vulnerability
  c_df$ses_vuln <- c_df$fw_status * c_df$ses_sens
  
  # Identify class breaks in vulnerability using Head/Tail scheme
  htb <- ht_breaks(x = c_df$ses_vuln, tsh = 0.8)
  
  
  # Generate maps using these class breaks ----
  
  # Adaptive capacity and ecological sensitivity per basin
  ac.f <- FeatureAreaAverage(FT.id = dis, 
                             RawDS = 1-ac15,
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
    ses.f <- 1 - (1 - ac.f)*(1-es.f) # fuzzy sum
  }
  if (Combs$agg[z] == 'Am') {
    ses.f <- (ac.f+es.f)/2 # mean 
  }
  
  # Classify each grid cell using Head/Tail breaks into transitional basins +, or not
  rclmat <- c(0, htb[1],        0, 
              htb[1], 1,        1) %>% 
    matrix(ncol = 3, byrow = T)
  transitional_and_hotspot_basins <- reclassify(ses.f*fw_status, rclmat, include.lowest = T)
  
  # Write classified raster for future analysis
  writeRaster(transitional_and_hotspot_basins, 
              paste0(here('Data/Sensitivity'), "/", 
                     Combs$disL[z], "_", 
                     Combs$agg[z], "_", 
                     Combs$useC[z], "_", 
                     Combs$sflw[z], "_transitional", sep = ""), 
              format = 'GTiff', overwrite = T)
  
  # Classify each grid cell using Head/Tail breaks into hotspot basin or not
  rclmat <- c(0, htb[2],        0, 
              htb[2], 1,        1) %>% 
    matrix(ncol = 3, byrow = T)
  hotspot_basins <- reclassify(ses.f*fw_status, rclmat, include.lowest = T)
  
  # Write classified raster for future analysis
  writeRaster(hotspot_basins, 
              paste0(here('Data/Sensitivity'), "/", 
                     Combs$disL[z], "_", 
                     Combs$agg[z], "_", 
                     Combs$useC[z], "_", 
                     Combs$sflw[z], "_hotspot", sep = ""),
              format = 'GTiff', overwrite = T)
  
  print(z)
  
} # End of loop