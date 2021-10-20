################################################################################
# @Manuscript - "Hotspots of social and ecological impacts from freshwater stress and storage loss" (Huggins et al.) 
# @Description - Derive freshwater stress and compare to freshwater storage trends per basin. Plot results. 
################################################################################

# Load 'here' library for easy path management.  
library(here)

# Import all setup and user-defined functions in R/setup and R/udfs folders
invisible(sapply(paste0(here("R/setup"), "/", list.files(here("R/setup"))), source)) 
invisible(sapply(paste0(here("R/udfs"), "/", list.files(here("R/udfs"))), source))

# Import data ----
wuse <- raster(here('Data/Dimensions/TotalWithdrawls_2010.tif'))
roff <- raster(here('Data/GSCD_Qmean_0d5.tif'))
tws  <- raster(here('Data/Rodell_GRACE_TWSt_raw.tif')) * 10 #cm to mm
dis  <- raster(here('Data/HyBas4_cleaned.tif')) # discretization scheme

# Calculate area-weighted basin average of all datasets
wuse <- FeatureAreaAverage(FT.id = dis, 
                           RawDS = wuse,
                           AreaDS = WGS84_areaRaster(0.5), 
                           operation = 'mean',
                           varnam = 'cons')

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
                       
# Reclassify freshwater stress, with class breaks at 10%, 40%
rclmat <- c(-Inf, 0.1, 0, 
            0.1, 0.4, 10,
            0.4, Inf, 20) %>% matrix(ncol = 3, byrow = T)
fwstrs <- reclassify(wuse/roff, rclmat) # Note reclassification is of ratio

# Reclassify TWS trends, with class breaks at -3, 3 mm/yr
rclmat <- c(-Inf, -3, 0,
            -3, 3, 1,
            3, Inf, 2) %>% matrix(ncol = 3, byrow = T)
twscls <- reclassify(tws, rclmat)

# Combine FW stress and TWS classes for bivariate plotting
strscomb <- fwstrs + twscls # these can be added, as ID values are at different magnitudes

# Plot three figures: (1) FW stress, (2) TWS trends, (3) Bivariate ----

# Import vector boundary files
disborders <- sf::read_sf(here('Data/hybas_l4_0d5.shp'))
coastlines <- sf::read_sf(here('Data/hybas_l4_coastlines.shp'))

# Reproject data for plotting in tmap
fwstrs.plt <- tmap_clipproj(wuse/roff) 
tws.plt <- tmap_clipproj(tws)
strscomb.plt <- tmap_clipproj(strscomb)

# Map of FW stress per basin
tm <- tm_shape(fwstrs.plt, projection = "+proj=robin") +
  tm_raster(style = "cont", palette = scico(20, palette = "bilbao",
                                            direction = 1), 
            midpoint = 0.2, breaks = c(0, 0.4), alpha = 1.00) +
  tm_shape(disborders) +
  tm_borders(lwd = 0.1, col = "grey40") +
  tm_shape(coastlines) +
  tm_borders(lwd = 0.7, col = "black") +
  tm_layout(legend.show = F, earth.boundary = c(-179, -60, 179, 88),
            earth.boundary.color = "white", space.color = "white",
            legend.frame = F, frame = F,
            outer.margins = c(-0, -0.09, -0, -0.04)) # B, L, T, R
tm    
tmap_save(tm, "C:/Users/xande/Desktop/jt_prep/1_stress_0_0d2_0d4.svg", units = "in")

# Map of TWS trend per basin
tm <- tm_shape(tws.plt, projection = "+proj=robin") +
  tm_raster(style = "cont", palette = scico(20, palette = "vikO",
                                            direction = -1)[3:18], 
            midpoint = 0, breaks = c(-20, 20), alpha = 1.00) +
  tm_shape(disborders) +
  tm_borders(lwd = 0.1, col = "grey40") +
  tm_shape(coastlines) +
  tm_borders(lwd = 0.7, col = "black") +
  tm_layout(legend.show = F, earth.boundary = c(-179, -60, 179, 88),
            earth.boundary.color = "white", space.color = "white",
            legend.frame = F, frame = F,
            outer.margins = c(-0, -0.09, -0, -0.04)) # B, L, T, R
tm
tmap_save(tm, "C:/Users/xande/Desktop/jt_prep/1_tws_dn20_0_d20.svg", units = "in")

# Map of bivariate combination of TWS and FW stress

# Custom color palette
bivar.pal <- c("#EEDAD9", "#E6E6E6", "#DEE8EB",
               "#BF9290", "#BFBFBF", "#9DB2B9",
               "#7E302E", "#7F7F7F", "#456A76")

tm <- tm_shape(strscomb.plt, projection = "+proj=robin") +
  tm_raster(style = "cat", palette = bivar.pal, alpha = 1.00) +
  tm_shape(disborders) +
  tm_borders(lwd = 0.1, col = "grey40") +
  tm_shape(coastlines) +
  tm_borders(lwd = 0.7, col = "black") +
  tm_layout(legend.show = F, earth.boundary = c(-179, -60, 179, 88),
            earth.boundary.color = "white", space.color = "white",
            legend.frame = F, frame = F,
            outer.margins = c(-0, -0.09, -0, -0.04)) # B, L, T, R
tm
tmap_save(tm, "C:/Users/xande/Desktop/jt_prep/1_bivar_0d10emph.svg", units = "in")


# Create bubble scatterplots for social-ecological activity distributions ----

# Re-import flux data as previously overwritten
wuse <- raster(here('Data/Dimensions/TotalWithdrawls_2010.tif'))
roff <- raster(here('Data/GSCD_Qmean_0d5.tif'))
tws  <- raster(here('Data/Rodell_GRACE_TWSt_raw.tif')) * 10 # cm to mm

# Import social-ecological dimension data
pop <- raster(here("Data/Dimensions/gpw_v4_population_count_adjusted_to_2015_unwpp_country_totals_rev11_2015_30_min.tif"))
cal <- raster(here("Data/Dimensions/SPAM_foodcrops_kcal_0d5.tif"))
gdp <- raster(here("Data/Dimensions/GDP_2015_0d5.tif"))
amphsr <- raster(here('Data/all_amph_0d5.tif'))
ramsar <- raster(here('Data/ramsar_count.tif'))

# Clean amphibian species richness of minor errors introduced through interpolation
amphsr[amphsr < 0 | amphsr > 129] <- 0

# Create dataframe of all dimensions
c_df <- raster::stack(dis, wuse, roff, 
                      tws, pop, cal, 
                      gdp, amphsr, ramsar, 
                      WGS84_areaRaster(0.5)) %>% as.data.frame() %>%
  set_colnames(c("dis", 'wuse', 'roff', 
                 'tws', 'pop', 'cal', 
                 'gdp', 'amphsr', 'ramsar', 'area'))
c_df <- c_df[complete.cases(c_df$dis),] # remove rows with NA discretization

# Calculate summary statistics per basin
c_df <- c_df %>% 
  group_by(dis) %>%
  summarize(
    tws  = weighted.mean(x = tws, w = area, na.rm = T),
    wuse = weighted.mean(x = wuse, w = area, na.rm = T),
    roff = weighted.mean(x = roff, w = area, na.rm = T),
    popE = sum(pop, na.rm = T),
    calE = sum(cal, na.rm = T),
    gdpE = sum(gdp, na.rm = T),
    amph = weighted.mean(x = amphsr, w = area, na.rm = T),
    ramsar = sum(ramsar, na.rm = T)
  )

# Calculate freshwater stress per basin
c_df$fwstrs <- rep(NA, nrow(c_df))
for (i in 1:nrow(c_df)) { c_df$fwstrs[i] <- c_df$wuse[i]/c_df$roff[i] }

# Remove basins that do not have adaptive capacity data 
disclude <- c(1130, 1520, 1840, 2170, 2360, 2450, 2530, 3550, 4240,
              4370, 5320, 5330, 5370, 5430, 5520, 5530, 5540, 5740,
              6540, 6730, 7220, 7610, 7630, 7650, 7710, 7740)

for (i in 1:length(disclude)) {
  c_df <- c_df %>% filter(dis != disclude[i])
}


# Clip TWS trends and FW stess values for plotting purposes 
for (i in 1:nrow(c_df)) {
  c_df$fwstrs[i] <- min(c_df$fwstrs[i], 1)
  c_df$tws[i] <- max(min(c_df$tws[i], 22), -22)
}

# Loop through plotting the four dimensions 
pltdims <- c('popE', 'calE', 'gdpE', 'amph')
aa = 0.4 # Transparency 
c_df$StressLog10 <- log10(100*c_df$fwstrs) # Convert stress to log10 scale 
c_df$StressLog10 <- ifelse(c_df$StressLog10 < 0, 0, c_df$StressLog10) # Handle small values

c_df$ramsar[c_df$ramsar == 0] <- NA # set 0s to NA so they don't appear when plotting
c_df$amph[c_df$amph == 0] <- NA # set 0s to NA so they don't appear when plotting

for (i in 1:length(pltdims)) {
  
  mp <- ggplot(data = c_df, aes_string(x = 'tws', y = 'StressLog10', size = pltdims[i] )) +
    annotate("rect", xmin=-Inf, xmax=-3, ymin=1, ymax=2, fill="#7D3232", alpha=aa) +  #upleft
    annotate("rect", xmin=-Inf, xmax=-3, ymin=0, ymax=1, fill="#7D6464", alpha=aa) +  #downleft
    annotate("rect", xmin=3, xmax=Inf, ymin=1,   ymax=2, fill="#19324B", alpha=aa) +  #upright
    annotate("rect", xmin=3, xmax=Inf, ymin=0,   ymax=1, fill="#4B647D", alpha=aa) +  #downright
    annotate("rect", xmin=-3, xmax=3, ymin=1,    ymax=2, fill="#E6E6E6", alpha=aa) +  #centre  
    geom_point(alpha= 0.5, shape = 21, color = 'black', stroke = 2) +
    {if(i != 4)scale_size(range = c(.1, 24))} +
    {if(i == 4)scale_size(range = c(.1, 12))} +
    {if(i == 4)
      geom_point(data = c_df, aes_string(x = 'tws', y = 'StressLog10', size = 'ramsar'),
                 alpha = 0.5, shape = 21, color = "#029E97", stroke = 2) } +
    xlab('') + ylab('')+
    coord_cartesian(xlim = c(-22, 22), ylim = c(0, 2), expand = c(0, 0), clip = "off") +
    scale_y_continuous(breaks = seq(0, 2, 1), 
                       limits = c(0, 2)) + 
    theme1 + theme(axis.ticks.x = element_line(size = 1),
                   axis.text = element_blank(),
                   plot.margin=unit(c(1,0.8,0.4,0.1),"cm"))  
  mp
  
  
  assign(paste0("dim_", pltdims[i], sep = ""), mp)
  
  ggsave(plot = last_plot(), 
         paste0("C:/Users/xande/Desktop/jt_prep/1_", pltdims[i], ".pdf", sep = ""),
         dpi = 500, width = 6, height = 4, units = "in")
}

# Calculate summary statistics ----

# Population 
sum(pop[strscomb == 0], na.rm =T)/1e9
sum(pop[strscomb == 10 | strscomb == 20], na.rm =T)/1e9

sum(pop[strscomb == 2], na.rm =T)/1e9
sum(pop[strscomb == 12 | strscomb == 22], na.rm =T)/1e9

# Verify with custom functions 
StressDry(c_df, 'popE', 1e9) #y
UnstressDry(c_df, 'popE', 1e9) #y
StressWet(c_df, 'popE', 1e9) #y
UnstressWet(c_df, 'popE', 1e9) #y

# Crop calories
StressDry(c_df, 'calE', 1e15)
UnstressDry(c_df, 'calE', 1e15)
StressWet(c_df, 'calE', 1e15)
UnstressWet(c_df, 'calE', 1e15)

# GDP
StressDry(c_df, 'gdpE', 1e12)
UnstressDry(c_df, 'gdpE', 1e12)
StressWet(c_df, 'gdpE', 1e12)
UnstressWet(c_df, 'gdpE', 1e12)

# Ramsar wetlands
StressDry(c_df, 'ramsar', 1)
UnstressDry(c_df, 'ramsar', 1)
StressWet(c_df, 'ramsar', 1)
UnstressWet(c_df, 'ramsar', 1)

# Amphibian species richness 
amphsr[is.na(strscomb)] <- NA
temp_sr <- amphsr
temp_sr[strscomb != 0] <- NA
weighted.mean(x = temp_sr, w = WGS84_areaRaster(0.5), na.rm = T)

temp_sr <- amphsr
temp_sr[strscomb != 10 & strscomb != 20] <- NA
weighted.mean(x = temp_sr, w = WGS84_areaRaster(0.5), na.rm = T)

temp_sr <- amphsr
temp_sr[strscomb != 2] <- NA
weighted.mean(x = temp_sr, w = WGS84_areaRaster(0.5), na.rm = T)

temp_sr <- amphsr
temp_sr[strscomb != 12 & strscomb != 22] <- NA
weighted.mean(x = temp_sr, w = WGS84_areaRaster(0.5), na.rm = T)

# Basin counts
c_df %>% filter(tws < -3 & fwstrs >= 0.4) %>% nrow()
c_df %>% filter(tws < -3 & fwstrs >= 0.1 & fwstrs < 0.4) %>% nrow()
c_df %>% filter(tws < -3 & fwstrs < 0.1) %>% nrow()

c_df %>% filter(tws >= -3 & tws <= 3 & fwstrs >= 0.4) %>% nrow()
c_df %>% filter(tws >= -3 & tws <= 3 &  fwstrs >= 0.1 & fwstrs < 0.4) %>% nrow()
c_df %>% filter(tws >= -3 & tws <= 3 &  fwstrs < 0.1) %>% nrow()

c_df %>% filter(tws > 3 & fwstrs >= 0.4) %>% nrow()
c_df %>% filter(tws > 3 & fwstrs >= 0.1 & fwstrs < 0.4) %>% nrow()
c_df %>% filter(tws >3 & fwstrs < 0.1) %>% nrow()