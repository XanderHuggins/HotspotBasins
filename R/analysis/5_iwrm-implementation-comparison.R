################################################################################
# @Manuscript - "Hotspots of social and ecological impacts from freshwater stress and storage loss" (Huggins et al.) 
# @Description - Compare vulnerability results with implementation levels of IWRM.
################################################################################

# Load 'here' library for easy path management.  
library(here)

# Import all setup and user-defined functions in R/setup and R/udfs folders
invisible(sapply(paste0(here("R/setup"), "/", list.files(here("R/setup"))), source)) 
invisible(sapply(paste0(here("R/udfs"), "/", list.files(here("R/udfs"))), source))

# Import data ----
cind <- raster(here('Data/fwss_ind_comb.tif'))
ac15 <- raster(here::here("Data", "AdaptiveCap2015.tif"))
efn <- raster(here("Data/Dimensions/EFN_0d5.tif"))
vsi <- raster(here("Data/Dimensions/VSI_aetw_0d5.tif"))
dis  <- raster(here('Data/HyBas4_cleaned.tif'))
SDG651_17 <- readr::read_csv(here('Data/SDG651_2017_baseline_cleaned.csv'))
SDG651_20 <- readr::read_csv(here('Data/SDG651_2020_stripped.csv'))
faoGAUL <- read_sf(here('Data/go_countryborders_faogaul2014Polygon.shp'))

# Create ID raster based on FAO GAUL ----
faoGAUL$uid <- seq(1, nrow(faoGAUL), 1)

# Raterize at fine resolution
GAULras <- fasterize(sf = faoGAUL, raster = WGS84_areaRaster(0.05), field = 'uid')

# Aggregate to operating resolution using modal value
GAULras <- aggregate(GAULras, fact = 10, expand = FALSE, fun = modal, na.rm = T)

# Convert SDG 6.5.1: IWRM implementation scores to raster, by joining w/ GAUL ----
SDGgaul_17 <- merge(x = faoGAUL, y = SDG651_17, by.x = 'ccode', by.y = 'ISOcode')
SDGgaul_20 <- merge(x = faoGAUL, y = SDG651_20, by.x = 'ccode', by.y = 'ISO3')

# Rasterize at fine resolution
SDGfinal_17 <- fasterize(sf = SDGgaul_17, raster = WGS84_areaRaster(0.05), field = 'Final', fun = 'min')

SDGfinal_20 <- fasterize(sf = SDGgaul_20, raster = WGS84_areaRaster(0.05), field = 'SDG651_2020', fun = 'min')

# Aggregate to operating resolution using mean value (ignore area-differences among contributing cells)
SDGfinal_17 <- aggregate(SDGfinal_17, fact = 10, expand = FALSE, fun = mean, na.rm = T)
SDGfinal_20 <- aggregate(SDGfinal_20, fact = 10, expand = FALSE, fun = mean, na.rm = T)

# Replace NAs in 2020 with values in 2017 if they exist
SDGfinal_latest <- SDGfinal_20
SDGfinal_latest[is.na(SDGfinal_20)] <- SDGfinal_17[is.na(SDGfinal_20)]

# Reproduce vulnerability inputs from script 4_... ----

# Social ecological sensitivity inputs, see comments in script 4_...
dis.ext <- dis
dis.ext[dis.ext >= 0] <- 1
ef.ptl <- RasterAreaPercentiles(RasterToClassify = efn,
                                WeightRaster = WGS84_areaRaster(0.5),
                                MaskRaster = dis.ext,
                                clipToExtent = "clip",
                                CRS.set = crs(dis.ext),
                                ext.set = extent(dis.ext))
ef.ptl[ef.ptl > 0] <- (1 - ef.ptl[ef.ptl > 0]) + 0.01 

vsi[is.na(vsi) | dis.ext != 1] <- NA
vsi.ptl <- RasterAreaPercentiles(RasterToClassify = vsi,
                                 WeightRaster = WGS84_areaRaster(0.5),
                                 MaskRaster = dis.ext,
                                 clipToExtent = "clip",
                                 CRS.set = crs(dis.ext),
                                 ext.set = extent(dis.ext))

# Inversion of adaptive capacity
ac <- 1 - ac15


# Create IWRM  mask to remove missing regions from analysis -----
Mask <- raster(SDGfinal_latest)
Mask[] <- 0
Mask[SDGfinal_latest > 0] <- 1


# Calculate basin average values of vulnerability inputs and IWRM implementation levels ----
c_df <- raster::stack(dis, cind, ef.ptl, vsi.ptl, ac, SDGfinal_latest, Mask, GAULras,
                      WGS84_areaRaster(0.5)) %>% as.data.frame() %>% 
  set_colnames(c('dis', 'cind', 'efn',  'vsi', 'ac', 'sdg651', 'mask', 'cc', 'area'))
c_df <- c_df[complete.cases(c_df$dis),]

# Remove basins that do not have adaptive capacity data 
disclude <- c(1130, 1520, 1840, 2170, 2360, 2450, 2530, 3550, 4240,
              4370, 5320, 5330, 5370, 5430, 5520, 5530, 5540, 5740,
              6540, 6730, 7220, 7610, 7630, 7650, 7710, 7740)

for (i in 1:length(disclude)) {
  c_df <- c_df %>% filter(dis != disclude[i])
}

c_df <- c_df %>% 
  group_by(dis) %>% 
  summarise(
    cind = weighted.quantile(x = cind, w = area, probs = 0.5, na.rm = T),
    efn = weighted.mean(x = efn, w = area, na.rm = T),
    vsi = weighted.mean(x = vsi, w = area, na.rm = T),
    ac = weighted.mean(x = ac, w = area, na.rm = T),
    sdg = weighted.mean(x = sdg651, w = area, na.rm = T),
    mask = weighted.mean(x = mask, w = area, na.rm = T),
    ccode = getmode(v = cc), # Identify most common country per basin
    ccode2 = getmode2(v = cc), # Identify second most common country per basin    
    tsbd = length(unique(cc)) # Identify number of countries per basin
  )

# Derive ecological sensitivity
c_df$ecosens <- (c_df$vsi + c_df$efn)/2
c_df$ecosens <- c_df$ecosens/max(c_df$ecosens, na.rm = T)

# Derive social-ecological sensitivity from fuzzy sum of eco. sens. and inv. ac. 
c_df$overallsens <- 1 - (1-c_df$ecosens)*(1-c_df$ac) # this is fuzzy sum

# Social-ecological vulnerabiltiy as product of sensitivity and basin freshwater status
c_df$overallprod <- c_df$cind * c_df$overallsens

# Identify vulnerability class breaks using Head/Tail breaks
htb_o <- ht_breaks2.0(x = c_df$overallprod, tsh = 0.8)

# Classify basins into hotspots 
c_df$ovhot <- rep(NA, nrow(c_df))
for (i in 1:nrow(c_df)) {
  c_df$ovhot[i] <- ifelse(c_df$overallprod[i] >= htb_o[2], 2, 1)
}

# Merge country codes with GAUL ID raster values, for labeling
gaulid <- data.frame(faoGAUL$ccode, faoGAUL$uid)
c_df <- merge(x = c_df, y = gaulid, by.x = 'ccode', by.y = 'faoGAUL.uid')

c_df$code2iso <- rep(NA, nrow(c_df)) 
for (i in 1:nrow(c_df)) {
  mode2 <- c_df$ccode2[i]
  ccodeiso <- gaulid %>% filter(faoGAUL.uid == mode2) %>% pull(faoGAUL.ccode)
  if (length(ccodeiso) == 0) {
    c_df$code2iso[i] <- NA  
  } else {
    c_df$code2iso[i] <- ccodeiso    
  }
}

c_df$label <- paste0(c_df$faoGAUL.ccode, ", ", c_df$code2iso, sep = "")

# Plot comparison of IWRM implementation levels and vulnerability results ----
aa = 0.4
c_df$rn <- seq(1, nrow(c_df), 1) # Create row number for periodic labeling
mp <- ggplot() +
  # background colors for IWRM class
  annotate("rect", xmin=0,  xmax=1,   ymin=0,  ymax=10,  fill="#8E292F", alpha=aa) + 
  annotate("rect", xmin=0,  xmax=1,   ymin=10, ymax=30,  fill="#D57F31", alpha=aa) + 
  annotate("rect", xmin=0,  xmax=1,   ymin=30, ymax=50,  fill="#E8BA24", alpha=aa) + 
  annotate("rect", xmin=0,  xmax=1,   ymin=50, ymax=70,  fill="#81BF81", alpha=aa) + 
  annotate("rect", xmin=0,  xmax=1,   ymin=70, ymax=90,  fill="#6CBDEF", alpha=aa) + 
  annotate("rect", xmin=0,  xmax=1,   ymin=90, ymax=100, fill="#0D7EB6", alpha=aa) + 
  geom_hline(yintercept = 10, lwd = 1, linetype = 'dashed') +
  geom_hline(yintercept = 30, lwd = 1, linetype = 'dashed') +
  geom_hline(yintercept = 50, lwd = 1, linetype = 'dashed') +
  geom_hline(yintercept = 70, lwd = 1, linetype = 'dashed') +
  geom_hline(yintercept = 90, lwd = 1, linetype = 'dashed') +
  geom_vline(xintercept = htb_o[1], lwd = 1, linetype = 'dashed') +
  geom_vline(xintercept = htb_o[2], lwd = 1, linetype = 'dashed') +
  geom_vline(xintercept = htb_o[3], lwd = 1, linetype = 'dashed') +
  # add red outline to transboundary basins
  geom_point(data = subset(c_df, mask >= 0.5 & tsbd >= 2), aes(x = overallprod, y = sdg), 
             alpha= 1, shape = 21, color = 'firebrick3', size = 3, lwd = 4) +
  
  geom_point(data = subset(c_df, mask >= 0.5 & tsbd >= 3), aes(x = overallprod, y = sdg), 
             alpha= 1, shape = 21, color = 'firebrick3', size = 4, lwd = 4) +
  
  geom_point(data = subset(c_df, mask >= 0.5 & tsbd >= 4), aes(x = overallprod, y = sdg), 
             alpha= 1, shape = 21, color = 'firebrick3', size = 5, lwd = 4) +
  
  geom_point(data = subset(c_df, mask >= 0.5 & tsbd >= 5), aes(x = overallprod, y = sdg), 
             alpha= 1, shape = 21, color = 'firebrick3', size = 6, lwd = 4) +
  
  geom_point(data = subset(c_df, mask >= 0.5 & tsbd >= 6), aes(x = overallprod, y = sdg), 
             alpha= 1, shape = 21, color = 'firebrick3', size = 7, lwd = 4) +
  # plot all basins satisfying mask constraint 
  
  # semi-opaque for non-transboundary
  geom_point(data = subset(c_df, mask >= 0.5 & tsbd == 1), aes(x = overallprod, y = sdg), 
             alpha= 0.25, shape = 19, color = 'black', size = 2) +
  
  # full opacity for transboundary
  geom_point(data = subset(c_df, mask >= 0.5 & tsbd >= 2), aes(x = overallprod, y = sdg), 
             alpha= 1, shape = 19, color = 'black', size = 2) +
  
  # mask areas under first break
  annotate("rect", xmin=0,  xmax=htb_o[2], ymin=0,  ymax=100, fill="white", alpha=1.5*aa) +
  scale_y_continuous(breaks = seq(10, 90, 20), 
                     limits = c(0, 100)) + 
  theme1 + theme(axis.ticks.x = element_line(size = 1)) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 100), expand = c(0, 0), clip = "on")
mp

ggsave(plot = mp, 
       paste0("C:/Users/xande/Desktop/jt_prep/5_iwrm_scatter_update.pdf", sep = ""),
       bg = 'transparent', dpi = 500, width = 6, height = 4, units = "in")


# Analysis of relationships ----
# Calculate correlation between overall vulnerability and IWRM implementation
rsq <- function (x, y) cor(x, y) ^ 2 # R squared function

# Evaluate for all basins
s_df <- c_df[complete.cases(c_df$sdg),]
rsq(s_df$overallprod, s_df$sdg/100)

# Evaluate for transitional or hotspot basins
s_df <- s_df %>% filter(overallprod >= htb_o[1])
rsq(s_df$overallprod, s_df$sdg/100)

# Evaluate for hotspot basins
s_df <- s_df %>% filter(overallprod >= htb_o[2])
rsq(s_df$overallprod, s_df$sdg/100)

# Identify number of basins 
c_df %>% nrow() ## output = 1204; total number of basins
c_df %>% filter(ovhot == 2) %>% nrow() ## output = 172; number of basins in hotspots
c_df %>% filter(ovhot == 2 & tsbd >= 2) %>% nrow() ## output = 61; # that are transboundary
61/172

c_df %>% filter(ovhot == 2 & sdg < 50) %>% nrow() ## output = 172; hotspot basins w/ low IWRM
c_df %>% filter(ovhot == 2 & tsbd >= 2 & sdg < 50) %>% nrow() ## output = 30

hot_trsb_sdg <- c_df %>% filter(ovhot == 2 & tsbd >= 2) %>% pull(sdg) 
hot_nontrsn_sdg <- c_df %>% filter(ovhot == 2 & tsbd < 2) %>% pull(sdg) 

t.test(x = hot_nontrsn_sdg, y = hot_trsb_sdg, alternative = "g")



# Create corresponding map ----

# Map inputs of adaptive capacity and ecological sensitivity
ac.f <- FeatureAreaAverage(FT.id = dis, RawDS = ac, # Already inverted
                           AreaDS = WGS84_areaRaster(0.5), operation = 'mean',
                           varnam = 'ac')

es.f <- FeatureAreaAverage(FT.id = dis, RawDS = (ef.ptl+vsi.ptl)/2,
                           AreaDS = WGS84_areaRaster(0.5), operation = 'mean',
                           varnam = 'ecosens')
# Scale ecological sensitivity, as required
es.f <- es.f/max(es.f[], na.rm = T)

# Social-ecological vulnerability
OvHot <- (1 - (1-ac.f)*(1-es.f)) * cind

# Reclassify results according to vulnerability classes
rclmat <- c(0, htb_o[2], 0,  
            htb_o[2], 1, 1) %>% matrix(ncol = 3, byrow = T)
OvHot_e <- reclassify(OvHot, rclmat)

# Generate emphasis polygons around hotspots
OvHot_e[OvHot_e == 0] <- NA
emph <- raster::rasterToPolygons(OvHot_e, na.rm = T, dissolve = T)


# Classify IWRM implementation levels according to IWRM Data Portal legend 
rclmat <- c(0,  10, 1,  
            10, 30, 2,
            30, 50, 3,
            50, 70, 4,
            70, 90, 5,
            90, 100,6) %>% matrix(ncol = 3, byrow = T)
IWRMclass <- reclassify(SDGfinal_latest, rclmat)


# Plot map 

# Import vector boundary files
disborders <- sf::read_sf(here('Data/hybas_l4_0d5.shp'))
coastlines <- sf::read_sf(here('Data/hybas_l4_coastlines.shp'))

# Reproject for Robinson projection
plt.obj <- tmap_clipproj(IWRMclass)

# Grey-out no data regions
grey.zn <- raster(IWRMclass)
grey.zn[OvHot >= 0 & is.na(IWRMclass)] <- 1
grey.zn <- tmap_clipproj(grey.zn)

# Darken non-hotspot regions
fade.zn <- raster(IWRMclass)
fade.zn[OvHot <= htb_o[2]] <- 1
fade.zn <- tmap_clipproj(fade.zn)

cus.pal <- c( #"#8E292F", 
  "#D57F31", "#E8BA24", "#81BF81", "#6CBDEF", "#0D7EB6")

tm <- 
  tm_shape(plt.obj, projection = "+proj=robin") + 
  tm_raster(style = "cat", palette = cus.pal) +
  tm_shape(grey.zn) +    tm_raster(style = "cat", palette = "grey80") +
  tm_shape(fade.zn) +    tm_raster(style = "cat", palette = "black", alpha = 0.4) +
  tm_shape(disborders) + tm_borders(lwd = 0.1, col = "grey40") +
  tm_shape(coastlines) + tm_borders(lwd = 0.7, col = "black") +
  tm_shape(emph) +       tm_borders(lwd = 2, col = "black") +
  tm_layout(legend.show = F, earth.boundary = c(-179, -60, 179, 88),
            earth.boundary.color = "white", space.color = "white",
            legend.frame = F, frame = F,
            outer.margins = c(-0, -0.09, -0, -0.04)) # B, L, T, R
tm

tmap_save(tm, "C:/Users/xande/Desktop/jt_prep/5_iwrm_update.svg", units = "in")