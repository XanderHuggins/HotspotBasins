# create plotting function

tmapsv <- function(var.to.plot, scicopalette, dir, lims, MP, svloc){
  disborders <- sf::read_sf(here('Data/hybas_l4_0d5.shp'))
  coastlines <- sf::read_sf(here('Data/hybas_l4_coastlines.shp'))
  
  # var.to.plot[var.to.plot == 0] <- NA
  var.plt <- tmap_clipproj(var.to.plot)
  
  # map of fw stress per fpu ----
  tm <- tm_shape(var.plt, projection = "+proj=robin") +
    tm_raster(style = "cont", palette = scico(20, palette = scicopalette,
                                              direction = dir), 
              breaks = lims, midpoint = MP, alpha = 1.00) +
    tm_shape(disborders) +
    tm_borders(lwd = 0.1, col = "grey40") +
    tm_shape(coastlines) +
    tm_borders(lwd = 0.7, col = "black") +
    tm_layout(legend.show = F, earth.boundary = c(-179, -60, 179, 88),
              earth.boundary.color = "white", space.color = "white",
              legend.frame = F, frame = F,
              outer.margins = c(-0, -0.09, -0, -0.04)) # B, L, T, R
  tm    
  tmap_save(tm, svloc, units = "in")
}



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

# Adaptive capacity - Varis et al. raw
tmapsv(var.to.plot = AdaptiveCap_2015, scicopalette = 'batlow', dir = -1, lims = c(0, 1), MP = 0.5,
       svloc = "C:/Users/xande/Desktop/ac_varis_raw.svg")

ac.f <- FeatureAreaAverage(FT.id = raster(here('Data/HyBas4_cleaned.tif')), 
                           RawDS = raster(here::here("Data", "AdaptiveCap2015.tif")),
                           AreaDS = WGS84_areaRaster(0.5), 
                           operation = 'mean',
                           varnam = 'ac')

tmapsv(var.to.plot = ac.f, scicopalette = 'batlow', dir = -1, lims = c(0, 1), MP = 0.5,
       svloc = "C:/Users/xande/Desktop/ac_basin_average.svg")




# gridded water withdrawals
wuse <- raster(here('Data/Dimensions/TotalWithdrawls_2010.tif'))
wuse[wuse <= 1] <- NA
tmapsv(var.to.plot = log10(wuse), scicopalette = 'lapaz', dir = -1, lims = c(0.1, 2),MP = 1, 
       svloc = "C:/Users/xande/Desktop/SImaps/gridwith.svg")


# basin average water withdrawals
wuse <- raster(here('Data/Dimensions/TotalWithdrawls_2010.tif'))
dis  <- raster(here('Data/HyBas4_cleaned.tif'))
wuse <- FeatureAreaAverage(FT.id = dis, RawDS = wuse,
                           AreaDS = WGS84_areaRaster(0.5), operation = 'mean',
                           varnam = 'cons')
wuse[wuse <= 1] <- NA
tmapsv(var.to.plot = log10(wuse), scicopalette = 'lapaz', dir = -1, lims = c(0.1, 2), MP = 1,
       svloc = "C:/Users/xande/Desktop/SImaps/baswith.svg")




# gridded tws trends
tws  <- raster(here('Data/Rodell_GRACE_TWSt_raw.tif')) * 10 #cm to mm

tmapsv(var.to.plot = tws, scicopalette = 'vikO', dir = -1, lims = c(-20, 20), MP = 0,
       svloc = "C:/Users/xande/Desktop/SImaps/gridtws.svg")


# basin average tws 
tws  <- raster(here('Data/Rodell_GRACE_TWSt_raw.tif')) * 10 #cm to mm
dis  <- raster(here('Data/HyBas4_cleaned.tif'))
tws <- FeatureAreaAverage(FT.id = dis, RawDS = tws,
                           AreaDS = WGS84_areaRaster(0.5), operation = 'mean',
                           varnam = 'tws')
tmapsv(var.to.plot = tws, scicopalette = 'vikO', dir = -1, lims = c(-20, 20), MP = 0,
       svloc = "C:/Users/xande/Desktop/SImaps/bastws.svg")




# gridded runoff
roff <- raster(here('Data/GSCD_Qmean_0d5.tif'))
tmapsv(var.to.plot = log10(roff), scicopalette = 'davos', dir = -1, lims = c(1, 3), MP = 2, 
       svloc = "C:/Users/xande/Desktop/SImaps/gridrunoff.svg")


# basin average runoff
roff <- raster(here('Data/GSCD_Qmean_0d5.tif'))
dis  <- raster(here('Data/HyBas4_cleaned.tif'))
roff <- FeatureAreaAverage(FT.id = dis, RawDS = roff,
                          AreaDS = WGS84_areaRaster(0.5), operation = 'mean',
                          varnam = 'roff')
tmapsv(var.to.plot = log10(roff), scicopalette = 'davos', dir = -1, lims = c(1, 3), MP = 2,
       svloc = "C:/Users/xande/Desktop/SImaps/basinrunoff.svg")





# plot basin freshwater stress
wuse <- raster(here('Data/Dimensions/TotalWithdrawls_2010.tif'))
roff <- raster(here('Data/GSCD_Qmean_0d5.tif'))

dis  <- raster(here('Data/HyBas4_cleaned.tif'))
wuse <- FeatureAreaAverage(FT.id = dis, RawDS = wuse,
                           AreaDS = WGS84_areaRaster(0.5), operation = 'mean',
                           varnam = 'cons')
roff <- FeatureAreaAverage(FT.id = dis, RawDS = roff,
                           AreaDS = WGS84_areaRaster(0.5), operation = 'mean',
                           varnam = 'cons')

tmapsv(var.to.plot = wuse/roff, scicopalette = 'bilbao', dir = 1, lims = c(0, 0.4), MP = 0.2,
       svloc = "C:/Users/xande/Desktop/SImaps/basfwSTRSS.svg")




# plot basin tws trends per runoff
tws  <- raster(here('Data/Rodell_GRACE_TWSt_raw.tif')) * 10 #cm to mm
roff <- raster(here('Data/GSCD_Qmean_0d5.tif'))

dis  <- raster(here('Data/HyBas4_cleaned.tif'))
tws <- FeatureAreaAverage(FT.id = dis, RawDS = tws,
                           AreaDS = WGS84_areaRaster(0.5), operation = 'mean',
                           varnam = 'cons')
roff <- FeatureAreaAverage(FT.id = dis, RawDS = roff,
                           AreaDS = WGS84_areaRaster(0.5), operation = 'mean',
                           varnam = 'cons')

tmapsv(var.to.plot = (tws/(0.4*roff)), scicopalette = 'vikO', dir = -1, lims = c(-1, 1), MP = 0,
       svloc = "C:/Users/xande/Desktop/SImaps/bastwsperQ.svg")




# plot combined indicator
# calculate area average of all inputs
wuse <- raster(here('Data/Dimensions/TotalWithdrawls_2010.tif'))
wuse <- FeatureAreaAverage(FT.id = dis, RawDS = wuse,
                           AreaDS = WGS84_areaRaster(0.5), operation = 'mean',
                           varnam = 'wuse')

roff <- raster(here('Data/GSCD_Qmean_0d5.tif'))
roff <- FeatureAreaAverage(FT.id = dis, RawDS = roff,
                           AreaDS = WGS84_areaRaster(0.5), operation = 'mean', 
                           varnam = 'roff')

tws  <- raster(here('Data/Rodell_GRACE_TWSt_raw.tif')) * 10 #cm to mm
tws <- FeatureAreaAverage(FT.id = dis, RawDS = tws, 
                          AreaDS = WGS84_areaRaster(0.5), operation = 'mean', 
                          varnam = 'tws')

# freshwater stress indicator
fws_ind <- min((wuse/(0.4*roff)), 1)

# storage trend indicator
tws_ind <- max(min((tws/(0.4*roff)), 1), -1)*-1
plot(tws_ind)

# combined indicator
cmb_ind <- max(min(( (fws_ind + tws_ind)/2 ), 1), 0)

tmapsv(var.to.plot = cmb_ind, scicopalette = 'lajolla', dir = 1, lims = c(0, 1), MP = 0.5,
       svloc = "C:/Users/xande/Desktop/SImaps/combind.svg")




# efn percentile
efn <- raster(here("Data/Dimensions/EFN_0d5.tif"))
dis.ext <- dis
dis.ext[dis.ext >= 0] <- 1

# convert all three into area-weighted percentiles, based on *contributing area*
ef.ptl <- RasterAreaPercentiles(RasterToClassify = efn,
                                WeightRaster = WGS84_areaRaster(0.5),
                                MaskRaster = dis.ext,
                                clipToExtent = "clip",
                                CRS.set = crs(dis.ext),
                                ext.set = extent(dis.ext))
ef.ptl[ef.ptl > 0] <- (1 - ef.ptl[ef.ptl > 0]) + 0.01 # invert as high depths are less sensitive

crs(ef.ptl) <- crs(dis)

tmapsv(var.to.plot = ef.ptl, scicopalette = 'batlow', dir = -1, lims = c(0, 1), MP = 0.5, 
       svloc = "C:/Users/xande/Desktop/SImaps/grid_efn_ptl.svg")


crs(efn) <- crs(dis)
tmapsv(var.to.plot = efn, scicopalette = 'batlow', dir = 1, lims = c(0, 2), MP = 1, 
       svloc = "C:/Users/xande/Desktop/SImaps/grid_efn.svg")


# vsi percentile
vsi <- raster(here("Data/Dimensions/VSI_aetw_0d5.tif"))
dis.ext <- dis
dis.ext[dis.ext >= 0] <- 1

# convert all three into area-weighted percentiles, based on *contributing area*
vsi[is.na(vsi) | dis.ext != 1] <- NA
vsi.ptl <- RasterAreaPercentiles(RasterToClassify = vsi,
                                 WeightRaster = WGS84_areaRaster(0.5),
                                 MaskRaster = dis.ext,
                                 clipToExtent = "clip",
                                 CRS.set = crs(dis.ext),
                                 ext.set = extent(dis.ext))

crs(vsi.ptl) <- crs(dis)

tmapsv(var.to.plot = vsi.ptl, scicopalette = 'batlow', dir = -1, lims = c(0, 1), MP = 0.5, 
       svloc = "C:/Users/xande/Desktop/SImaps/grid_vsi_ptl.svg")


crs(vsi) <- crs(dis)
tmapsv(var.to.plot = vsi/max(vsi[], na.rm = T), scicopalette = 'batlow', dir = -1, lims = c(0, 1), MP = 0.5, 
       svloc = "C:/Users/xande/Desktop/SImaps/grid_vsi.svg")



# Basin average vsi 
vsi <- raster(here("Data/Dimensions/VSI_aetw_0d5.tif"))
dis.ext <- dis
dis.ext[dis.ext >= 0] <- 1

# convert all three into area-weighted percentiles, based on *contributing area*
vsi[is.na(vsi) | dis.ext != 1] <- NA
vsi.ptl <- RasterAreaPercentiles(RasterToClassify = vsi,
                                 WeightRaster = WGS84_areaRaster(0.5),
                                 MaskRaster = dis.ext,
                                 clipToExtent = "clip",
                                 CRS.set = crs(dis.ext),
                                 ext.set = extent(dis.ext))
crs(vsi.ptl) <- crs(dis)

vsi.ptl <- FeatureAreaAverage(FT.id = dis, RawDS = vsi.ptl,
                           AreaDS = WGS84_areaRaster(0.5), operation = 'mean',
                           varnam = 'vsi.ptl')

tmapsv(var.to.plot = vsi.ptl, scicopalette = 'batlow', dir = -1, lims = c(0, 1), MP = 0.5, 
       svloc = "C:/Users/xande/Desktop/SImaps/basin_vsi_ptl.svg")



# basin average efn sensitivity
efn <- raster(here("Data/Dimensions/EFN_0d5.tif"))
dis.ext <- dis
dis.ext[dis.ext >= 0] <- 1

# convert all three into area-weighted percentiles, based on *contributing area*
ef.ptl <- RasterAreaPercentiles(RasterToClassify = raster(here("Data/Dimensions/EFN_0d5.tif")),
                                WeightRaster = WGS84_areaRaster(0.5),
                                MaskRaster = dis.ext,
                                clipToExtent = "clip",
                                CRS.set = crs(dis.ext),
                                ext.set = extent(dis.ext))
ef.ptl[ef.ptl > 0] <- (1 - ef.ptl[ef.ptl > 0]) + 0.01 # invert as high depths are less sensitive

crs(ef.ptl) <- crs(dis)

ef.ptl <- FeatureAreaAverage(FT.id = dis, RawDS = ef.ptl,
                              AreaDS = WGS84_areaRaster(0.5), operation = 'mean',
                              varnam = 'ef.ptl')

tmapsv(var.to.plot = ef.ptl, scicopalette = 'batlow', dir = -1, lims = c(0, 1), MP = 0.5, 
       svloc = "C:/Users/xande/Desktop/SImaps/basin_efn_ptl.svg")


# combined sensitivity indicator
sensInd <- (vsi.ptl+ef.ptl)/2
sensInd <- sensInd/max(sensInd[], na.rm = T)
tmapsv(var.to.plot = sensInd, scicopalette = 'batlow', dir = -1, lims = c(0, 1), MP = 0.5, 
       svloc = "C:/Users/xande/Desktop/SImaps/basin_sensitivity_indicator.svg")




## SI Figure 7 - Social-ecological sensitivity ----

dis.ext <- raster(here('Data/HyBas4_cleaned.tif'))
dis.ext[dis.ext >= 0] <- 1

ef.ptl <- RasterAreaPercentiles(RasterToClassify = raster(here("Data/Dimensions/EFN_0d5.tif")),
                                WeightRaster = WGS84_areaRaster(0.5),
                                MaskRaster = dis.ext,
                                clipToExtent = "clip",
                                CRS.set = crs(dis.ext),
                                ext.set = extent(dis.ext))
ef.ptl[ef.ptl > 0] <- (1 - ef.ptl[ef.ptl > 0]) + 0.01 # Invert as greater depths are less sensitive

# Vegetation productivity sensitivity to soil moisture & shallow groundwater anomalies
vsi <- raster(here("Data/Dimensions/VSI_aetw_0d5.tif"))
vsi[is.na(vsi) | dis.ext != 1] <- NA
vsi.ptl <- RasterAreaPercentiles(RasterToClassify = vsi,
                                 WeightRaster = WGS84_areaRaster(0.5),
                                 MaskRaster = dis.ext,
                                 clipToExtent = "clip",
                                 CRS.set = crs(dis.ext),
                                 ext.set = extent(dis.ext))

ac.f <- FeatureAreaAverage(FT.id = raster(here('Data/HyBas4_cleaned.tif')), 
                           RawDS = 1 - raster(here::here("Data", "AdaptiveCap2015.tif")),
                           AreaDS = WGS84_areaRaster(0.5), 
                           operation = 'mean',
                           varnam = 'ac')

# Calculate ecological sensitivity per basin
es.f <- (ef.ptl+vsi.ptl)/2
es.f <- FeatureAreaAverage(FT.id = raster(here('Data/HyBas4_cleaned.tif')), 
                           RawDS = es.f,
                           AreaDS = WGS84_areaRaster(0.5), 
                           operation = 'mean',
                           varnam = 'ecosens')

# Normalize by maximum basin value
es.f <- es.f/max(es.f[], na.rm = T)

# Calculate social-ecological sensitivity through fuzzy sum of inputs
ov.f <- 1 - (1 - ac.f)*(1-es.f) 

crs(es.f) <- crs(WGS84_areaRaster(0.5))
# plot each 
tmapsv(var.to.plot = es.f, scicopalette = 'batlow', dir = -1, lims = c(0, 1), MP = 0.5,
       svloc = "C:/Users/xande/Desktop/eco_sens.svg")

tmapsv(var.to.plot = ac.f, scicopalette = 'batlow', dir = -1, lims = c(0, 1), MP = 0.5, 
       svloc = "C:/Users/xande/Desktop/soc_sens.svg")

tmapsv(var.to.plot = ov.f, scicopalette = 'batlow', dir = -1, lims = c(0, 1), MP = 0.5, 
       svloc = "C:/Users/xande/Desktop/SES_sens.svg")



# basin inverted adaptive capacity
ac15 <- raster(here::here("Data", "AdaptiveCap2015.tif"))
ac15 <- FeatureAreaAverage(FT.id = dis, RawDS = 1-ac15,
                             AreaDS = WGS84_areaRaster(0.5), operation = 'mean',
                             varnam = 'ac15')


tmapsv(var.to.plot = ac15, scicopalette = 'batlow', dir = -1, lims = c(0, 1), MP = 0.5, 
       svloc = "C:/Users/xande/Desktop/SImaps/ac_inv_basin.svg")


# plot overall social+ecological sensitivity
OvInd <- 1- (1-sensInd)*(1-ac15)

tmapsv(var.to.plot = OvInd, scicopalette = 'batlow', dir = -1, lims = c(0, 1), MP = 0.5, 
       svloc = "C:/Users/xande/Desktop/SImaps/overall_sens.svg")
