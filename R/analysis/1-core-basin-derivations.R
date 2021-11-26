# Name: 1-core-basin-derivations.R
# Description: Derive basin level results (see list below) for all core analyses of the study 
# Core basin results derived in this script:
# - Basin freshwater status
# - Social, ecological, and social-ecological sensitivity
# - Social, ecological, and social-ecological vulnerability
# - Hotspot basin classification  
# - Population, food crop production, gross domestic product, amphibian species richness, and Ramsar Sites within each basin

# Use 'here' package for easy path management
library(here)

# Import all setup and user-defined functions in R/setup and R/udfs folders
invisible(sapply(paste0(here("R/setup"), "/", list.files(here("R/setup"))), source)) 
invisible(sapply(paste0(here("R/udf"), "/", list.files(here("R/udf"))), source))

# Import data ----
# Geospatial areal units of analysis
dis  <- raster(here('Data/HyBas4_cleaned.tif')) # basin scheme

# Hydrological data
wuse <- raster(here('Data/Dimensions/TotalWithdrawls_2010.tif'))
roff <- raster(here('Data/GSCD_Qmean_0d5.tif'))
tws  <- raster(here('Data/Rodell_GRACE_TWSt_raw.tif')) * 10 #cm to mm

# Adaptive capacity data
ac15 <- raster(here::here("Data", "AdaptiveCap2015.tif"))

# Ecological sensitivity data
efn <- raster(here("Data/Dimensions/EFN_0d5.tif"))
vsi <- raster(here("Data/Dimensions/VSI_aetw_0d5.tif"))

# Social-ecological activity data
pop <- raster(here("Data/Dimensions/gpw_v4_population_count_adjusted_to_2015_unwpp_country_totals_rev11_2015_30_min.tif"))
cal <- raster(here("Data/Dimensions/SPAM_foodcrops_kcal_0d5.tif"))
gdp <- raster(here("Data/Dimensions/GDP_2015_0d5.tif"))
amphsr <- raster(here('Data/all_amph_0d5.tif'))
ramsar <- raster(here('Data/ramsar_count.tif'))
amphsr[amphsr < 0 | amphsr > 129] <- 0 # Remove artifacts from interpolation 
wcc <- raster(here::here("Data", "Dimensions", "WCC_since2000.tif"))

# Transform data ----
# Convert environmental flow sensitivity and vegetation sensitivity rasters into area-weighted percentiles # Establish constant extent for masking 
dis.ext <- dis
dis.ext[dis.ext >= 0] <- 1

ef.ptl <- RasterAreaPercentiles(RasterToClassify = efn,
                                WeightRaster = WGS84_areaRaster(0.5),
                                MaskRaster = dis.ext,
                                clipToExtent = "clip",
                                CRS.set = crs(dis.ext),
                                ext.set = extent(dis.ext))
ef.ptl[ef.ptl > 0] <- (1 - ef.ptl[ef.ptl > 0]) + 0.01 # Invert as greater depths are less sensitive

# Vegetation productivity sensitivity to soil moisture & shallow groundwater anomalies
vsi[is.na(vsi) | dis.ext != 1] <- NA
vsi.ptl <- RasterAreaPercentiles(RasterToClassify = vsi,
                                 WeightRaster = WGS84_areaRaster(0.5),
                                 MaskRaster = dis.ext,
                                 clipToExtent = "clip",
                                 CRS.set = crs(dis.ext),
                                 ext.set = extent(dis.ext))

# Stack all data and convert to dataframe for faster processing, and calculate summary statistics ----
c_df <- raster::stack(dis, # basins
                      wuse, roff, tws, # hydrological data 
                      ef.ptl, vsi.ptl, ac15, # Adaptive capacity and ecological sensitivity data
                      pop, cal, gdp, amphsr, ramsar,  # Social-ecological activity 
                      wcc, # World Conflict Chronology data
                      WGS84_areaRaster(0.5)) %>% 
  as.data.frame() %>% 
  set_colnames(c('id', 
                 'wuse', 'roff', 'tws',
                 'efn',  'vsi', 'ac', 
                 'pop', 'cal', 'gdp', 'amph','ramsar', 
                 'wcc', 
                 'area'))
c_df <- c_df[complete.cases(c_df$id),]

c_df <- c_df %>% 
  group_by(id) %>%
  summarize(
    tws  = weighted.mean(x = tws, w = area, na.rm = T),
    wuse = weighted.mean(x = wuse, w = area, na.rm = T),
    roff = weighted.mean(x = roff, w = area, na.rm = T),
    efn = weighted.mean(x = efn, w = area, na.rm = T),
    vsi = weighted.mean(x = vsi, w = area, na.rm = T),
    ac = weighted.mean(x = ac, w = area, na.rm = T),
    pop = sum(pop, na.rm = T)/1e9,
    cal = sum(cal, na.rm = T)/1e15,
    gdp = sum(gdp, na.rm = T)/1e12,
    amph = weighted.mean(x = amph, w = area, na.rm = T),
    ramsar = sum(ramsar, na.rm = T),
    wcc = sum(wcc, na.rm = T),
    area = sum(area, na.rm = T)
  )

# Calculate basin freshwater status ----

# Calculate freshwater stress per basin
c_df$fwstrs <- rep(NA)
for (i in 1:nrow(c_df)) { c_df$fwstrs[i] <- c_df$wuse[i]/c_df$roff[i] }

# Freshwater stress indicator
c_df$fws_ind <- rep(NA)
for (i in 1:nrow(c_df)) { c_df$fws_ind[i] <- min(c_df$wuse[i]/(0.4*c_df$roff[i]), 1) }

# Freshwater storage trend indicator
c_df$tws_ind <- rep(NA)
for (i in 1:nrow(c_df)) { c_df$tws_ind[i] <- max(min((c_df$tws[i]/(0.4*c_df$roff[i])), 1), -1) * -1 }

# Basin freshwater status
c_df$fw_status <- rep(NA)
for (i in 1:nrow(c_df)) { c_df$fw_status[i] <- max(min(( (c_df$fws_ind[i] + c_df$tws_ind[i])/2 ), 1), 0) }

# Handle basins with earthquake interference by setting status to the freshwater stress indicator

# List of basins identified with EQ interference
eq_intfr <- readr::read_csv(file = here::here("Data/basin_lists_1.csv"),
                           col_select = 'eq_interference',
                           col_names = T)

for (i in 1:nrow(c_df)) { 
  if (c_df$id[i] %in% eq_intfr$eq_interference) {
    c_df$fw_status[i] <- c_df$fws_ind[i]
    }
}

# Calculate ecological sensitivity indicator  ----
c_df$ecol_sens <- (c_df$vsi + c_df$efn)/2
c_df$ecol_sens <- c_df$ecol_sens/max(c_df$ecol_sens, na.rm = T)

# Show that inverted adaptive capacity is used as social sensitivity indicator ---- 
c_df$socl_sens <- 1 - c_df$ac

# Derive social-ecological sensitivity as fuzzy sum of ecological and social sensitivities ----
c_df$ses_sens <- 1 - ((1 - c_df$ecol_sens) * (1- c_df$socl_sens))

# Derive social, ecological, and social-ecological vulnerability as the product of sensitivity and the freshwater basin status ----
c_df$ecol_vuln <- c_df$ecol_sens * c_df$fw_status
c_df$socl_vuln <- c_df$socl_sens * c_df$fw_status
c_df$ses_vuln  <- c_df$ses_sens  * c_df$fw_status

# Identify vulnerability classes using the head/tails breaks method for each vulnerability 
htb_ecol <- ht_breaks(x = c_df$ecol_vuln, tsh = 0.8) # threshold does not effect outcomes, but ensures min. 3 breaks are identified for each 
htb_socl <- ht_breaks(x = c_df$socl_vuln, tsh = 0.8)
htb_ses <- ht_breaks(x = c_df$ses_vuln, tsh = 0.8)

# Using these breaks, classify each dimension into vulnerability classes
c_df$ecol_vclass <- rep(NA)
c_df$socl_vclass <- rep(NA)
c_df$ses_vclass  <- rep(NA)

# Ecological classes
c_df$ecol_vclass[c_df$ecol_vuln <  htb_ecol[1]] <- 0 # Low vulnerability
c_df$ecol_vclass[c_df$ecol_vuln >= htb_ecol[1]] <- 1 # Moderate vulnerability  
c_df$ecol_vclass[c_df$ecol_vuln >= htb_ecol[2]] <- 2 # High vulnerability 
c_df$ecol_vclass[c_df$ecol_vuln >= htb_ecol[3]] <- 3 # Very high vulnerability 

# Social classes
c_df$socl_vclass[c_df$socl_vuln <  htb_socl[1]] <- 0 # Low vulnerability
c_df$socl_vclass[c_df$socl_vuln >= htb_socl[1]] <- 1 # Moderate vulnerability  
c_df$socl_vclass[c_df$socl_vuln >= htb_socl[2]] <- 2 # High vulnerability 
c_df$socl_vclass[c_df$socl_vuln >= htb_socl[3]] <- 3 # Very high vulnerability 

c_df$ses_vclass[c_df$ses_vuln <  htb_ses[1]] <- 0 # Low vulnerability
c_df$ses_vclass[c_df$ses_vuln >= htb_ses[1]] <- 1 # Moderate vulnerability  
c_df$ses_vclass[c_df$ses_vuln >= htb_ses[2]] <- 2 # HOTSPOT BASIN - High vulnerability 
c_df$ses_vclass[c_df$ses_vuln >= htb_ses[3]] <- 3 # HOTSPOT BASIN - Very high vulnerability 
  

# Merge with basin shapefile and write for plotting and statistical analyses
basins <- sf::read_sf(here('Data/hybas_l4_0d5.shp'))
names(basins) <- c('pfaf_id', 'geometry')
basins <- merge(x = basins, y = c_df, by.x = 'pfaf_id', by.y = 'id', all.x = TRUE)

write_sf(obj = basins, dsn = here('R/plotting+stats/'), layer = 'Basin_data.shp', 
         driver = 'ESRI Shapefile')
