################################################################################
# @Manuscript - "Hotspots of social and ecological impacts from freshwater stress and storage loss" (Huggins et al.) 
# @Description - Compare basin freshwater status with social adaptive capacity. Plot results.
################################################################################

# Load 'here' library for easy path management.  
library(here)

# Import all setup and user-defined functions in R/setup and R/udfs folders
invisible(sapply(paste0(here("R/setup"), "/", list.files(here("R/setup"))), source)) 
invisible(sapply(paste0(here("R/udfs"), "/", list.files(here("R/udfs"))), source))

# Import data ----
cind <- raster(here('Data/fwss_ind_comb.tif'))
ac15 <- raster(here::here("Data", "AdaptiveCap2015.tif"))
ac15[is.na(cind)] <- NA
dis  <- raster(here('Data/HyBas4_cleaned.tif'))

# Calculate basin statistics before plotting as percentile thresholds are used ----
c_df <- raster::stack(dis, cind, ac15, WGS84_areaRaster(0.5)) %>% 
  as.data.frame() %>% set_colnames(c('dis', 'cind', 'ac', 'area'))
c_df <- c_df[complete.cases(c_df$dis),]

# Remove basins that do not have adaptive capacity 
disclude <- c(1130, 1520, 1840, 2170, 2360, 2450, 2530, 3550, 4240,
              4370, 5320, 5330, 5370, 5430, 5520, 5530, 5540, 5740,
              6540, 6730, 7220, 7610, 7630, 7650, 7710, 7740)

for (i in 1:length(disclude)) {
  c_df <- c_df %>% filter(dis != disclude[i])
}

s_df <- c_df %>% 
  group_by(dis) %>% 
  summarize(
    cind = weighted.quantile(x = cind, w = area, probs = 0.5, na.rm = T), # Median as basin freshwater status is already provided at the basin scale -- this ensures any unwanted edge effects do not change basin values
    ac   = weighted.mean(x = ac, w = area, na.rm = T)
  )

# Calculate percentile values
ptls <- data.frame(P = seq(1, 100, by = 1), 
                   ac  = rep(NA, 100),
                   ind = rep(NA, 100))

for (i in 1:nrow(ptls)) {
  ptls$ac[i]  <- quantile(x = s_df$ac,   probs = ptls$P[i]/100, na.rm = T)
  ptls$ind[i] <- quantile(x = s_df$cind, probs = ptls$P[i]/100, na.rm = T)
}

# Create raster of adaptive capacity per basin
ac15 <- FeatureAreaAverage(FT.id = dis, 
                           RawDS = ac15,
                           AreaDS = WGS84_areaRaster(0.5), 
                           operation = 'mean',
                           varnam = 'ac')

# Reclassify basin freshwater status and adaptive capacity ----
# Adaptive capacity
rclmat <- c(0, ptls$ac[33], 0,  
            ptls$ac[33], ptls$ac[67], 10, 
            ptls$ac[67], 1, 20) %>% matrix(ncol = 3, byrow = T)
ac_cls <- reclassify(ac15, rclmat)

# Basin freshwater status
rclmat <- c(0, ptls$ind[67], 0, # P67
            ptls$ind[67], ptls$ind[80], 1, #P80
            ptls$ind[80], 1, 2) %>% matrix(ncol = 3, byrow = T)
ci_cls <- reclassify(cind, rclmat)

# Combine reclassified rasters for bivariate plotting 
ses_comb <- ac_cls + ci_cls # these can be added, as ID values are at different magnitudes

writeRaster(ses_comb, here('Data/Fig2_combs.tif'), format = 'GTiff',
            overwrite = T)

c_df <- raster::stack(ses_comb, dis, WGS84_areaRaster(0.5)) %>% 
  as.data.frame() %>% set_colnames(c('bivar', 'dis', 'area'))
c_df <- c_df[complete.cases(c_df$dis),]

s_df <- c_df %>% 
  group_by(bivar) %>% 
  summarize(
    basin_count = length(unique(dis))
  )
s_df

# Plot two figures: (1) Adaptive capacity, (2) Bivariate ----

# Import vector boundary files
disborders <- sf::read_sf(here('Data/hybas_l4_0d5.shp'))
coastlines <- sf::read_sf(here('Data/hybas_l4_coastlines.shp'))

cus.pal <- c("#7F7F7F", "#19255A", "#7E1900",# bottom row of bivar legend, L to R
             "#BFBFBF", "#455186", "#9D5819",
             "#E6E6E6", "#7683B5", "#B99232")# upper row of bivar legend, L to R

# Create adaptive capacity plot
plt.obj <- tmap_clipproj(1-ac15) # reproject for plotting

tm <- tm_shape(plt.obj, projection = "+proj=robin") +
  tm_raster(style = "cont", palette = scico(20, palette = "batlow", direction = 1),
            breaks= c(0,1), midpoint = 0.5) +
  tm_shape(disborders) +
  tm_borders(lwd = 0.1, col = "grey40") +
  tm_shape(coastlines) +
  tm_borders(lwd = 0.7, col = "black") +
  tm_layout(legend.show = F, earth.boundary = c(-179, -60, 179, 88),
            earth.boundary.color = "white", space.color = "white",
            legend.frame = F, frame = F,
            outer.margins = c(-0, -0.09, -0, -0.04)) # B, L, T, R
tm
tmap_save(tm, 
          "C:/Users/xande/Desktop/jt_prep/2_ac.svg", units = "in")

# Create bivariate plot
plt.obj <- tmap_clipproj(ses_comb) # reproject for plotting

tm <- tm_shape(plt.obj, projection = "+proj=robin") +
  tm_raster(style = "cat", palette = cus.pal) +
  tm_shape(disborders) +
  tm_borders(lwd = 0.1, col = "grey40") +
  tm_shape(coastlines) +
  tm_borders(lwd = 0.7, col = "black") +
  tm_layout(legend.show = F, earth.boundary = c(-179, -60, 179, 88),
            earth.boundary.color = "white", space.color = "white",
            legend.frame = F, frame = F,
            outer.margins = c(-0, -0.09, -0, -0.04)) # B, L, T, R
tm
tmap_save(tm, 
          "C:/Users/xande/Desktop/jt_prep/2_ac_ind.svg", units = "in")


# Generate scatterplots for social-ecological activity distributions ----

# Re-import indicator data as previously overwritten
cind <- raster(here('Data/fwss_ind_comb.tif'))
ac15 <- raster(here::here("Data", "AdaptiveCap2015.tif")) # invert on import

# Import social-ecological dimension data
pop <- raster(here("Data/Dimensions/gpw_v4_population_count_adjusted_to_2015_unwpp_country_totals_rev11_2015_30_min.tif"))
cal <- raster(here("Data/Dimensions/SPAM_foodcrops_kcal_0d5.tif"))
gdp <- raster(here("Data/Dimensions/GDP_2015_0d5.tif"))
amphsr <- raster(here('Data/all_amph_0d5.tif'))
ramsar <- raster(here('Data/ramsar_count.tif'))

# Clean amphibian species richness of minor errors introduced through interpolation
amphsr[amphsr < 0 | amphsr > 129] <- 0

# Create dataframe of all dimensions
c_df <- raster::stack(dis, cind, ac15, pop, cal, gdp, amphsr, ramsar,
                   WGS84_areaRaster(0.5)) %>% as.data.frame() %>% 
  set_colnames(c("dis", 'cind', 'ac', 'pop', 'cal', 'gdp', 'amphsr', 'ramsar', 'area'))
c_df <- c_df[complete.cases(c_df$dis),] # remove NA dis rows

# Remove basins that do not have adaptive capacity data 
disclude <- c(1130, 1520, 1840, 2170, 2360, 2450, 2530, 3550, 4240,
              4370, 5320, 5330, 5370, 5430, 5520, 5530, 5540, 5740,
              6540, 6730, 7220, 7610, 7630, 7650, 7710, 7740)

for (i in 1:length(disclude)) {
  c_df <- c_df %>% filter(dis != disclude[i])
}

# Calculate summary statistics per basin
c_df <- c_df %>% 
  group_by(dis) %>%
  summarize(
    cind = weighted.quantile(x = cind, w = area, probs = 0.5, na.rm = T),
    ac   = weighted.mean(x = ac, w = area, na.rm = T),
    popE = sum(pop, na.rm = T),
    calE = sum(cal, na.rm = T),
    gdpE = sum(gdp, na.rm = T),
    amph = weighted.mean(x = amphsr, w = area, na.rm = T),
    ramsar = sum(ramsar, na.rm = T)
  )

# Loop through plotting the four dimensions
pltdims <- c('popE', 'calE', 'gdpE', 'amph')
aa = 0.4

c_df$ramsar[c_df$ramsar == 0] <- NA # set 0s to NA so they don't appear when plotting
c_df$amph[c_df$amph == 0] <- NA # set 0s to NA so they don't appear when plotting

for (i in 1:length(pltdims)) {
  
  mp <- ggplot(data = c_df, aes_string(x = 'cind', y = 'ac', size = pltdims[i] )) +
    annotate("rect", xmin=0,            xmax=ptls$ind[67], ymin=ptls$ac[67], ymax=1,           fill="#E6E6E6", alpha=aa) +  
    annotate("rect", xmin=0,            xmax=ptls$ind[67], ymin=ptls$ac[33], ymax=ptls$ac[67], fill="#BFBFBF", alpha=aa) +  
    annotate("rect", xmin=0,            xmax=ptls$ind[67], ymin=0,           ymax=ptls$ac[33], fill="#7F7F7F", alpha=aa) +     
    annotate("rect", xmin=ptls$ind[67], xmax=ptls$ind[80], ymin=ptls$ac[67], ymax=1,           fill="#7683B5", alpha=aa) +     
    annotate("rect", xmin=ptls$ind[67], xmax=ptls$ind[80], ymin=ptls$ac[33], ymax=ptls$ac[67], fill="#455186", alpha=aa) +
    annotate("rect", xmin=ptls$ind[67], xmax=ptls$ind[80], ymin=0,           ymax=ptls$ac[33], fill="#19255A", alpha=aa) +
    annotate("rect", xmin=ptls$ind[80], xmax=1,            ymin=ptls$ac[67], ymax=1,           fill="#B99232", alpha=aa) +  
    annotate("rect", xmin=ptls$ind[80], xmax=1,            ymin=ptls$ac[33], ymax=ptls$ac[67], fill="#9D5819", alpha=aa) +  
    annotate("rect", xmin=ptls$ind[80], xmax=1,            ymin=0,           ymax=ptls$ac[33], fill="#7E1900", alpha=aa) +
    geom_point(alpha= 0.5, shape = 21, color = 'black', stroke = 2) +
    {if(i != 4)scale_size(range = c(0.1, 24))} +
    {if(i == 4)scale_size(range = c(0.1, 12))} +
    {if(i == 4)
      geom_point(data = c_df, aes_string(x = 'cind', y = 'ac', size = 'ramsar'),
                 alpha = 0.5, shape = 21, color = "#029E97", stroke = 2) } +
    xlab('') + ylab('')+
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = c(0, 0), clip = "off") +
    scale_y_continuous(breaks = seq(0, 1, 0.2), 
                       limits = c(0, 1)) + 
    theme1 + theme(axis.ticks.x = element_line(size = 1),
                   axis.text = element_blank(),
                   plot.margin=unit(c(1,0.8,0.4,0.6),"cm"))   
  
  mp
  
  assign(paste0("dim_", pltdims[i], sep = ""), mp)
  
  ggsave(plot = last_plot(), 
         paste0("C:/Users/xande/Desktop/jt_prep/2_", pltdims[i], ".pdf", sep = ""),
         dpi = 500, width = 6, height = 4, units = "in")
}


# Calculate summary statistics ----

# Population counts and percentages
sum(pop[ses_comb == 1], na.rm =T)/1e9
sum(pop[ses_comb == 1], na.rm =T)/sum(pop[],na.rm = T)

sum(pop[ses_comb == 11 | ses_comb == 21], na.rm =T)/1e9
sum(pop[ses_comb == 11 | ses_comb == 21], na.rm =T)/sum(pop[],na.rm = T)

sum(pop[ses_comb == 2], na.rm =T)/1e9
sum(pop[ses_comb == 2], na.rm =T)/sum(pop[],na.rm = T)

sum(pop[ses_comb == 12 | ses_comb == 22], na.rm =T)/1e9
sum(pop[ses_comb == 12 | ses_comb == 22], na.rm =T)/sum(pop[],na.rm = T)

# Food crop counts and percentages
sum(cal[ses_comb == 1], na.rm =T)/1e15
sum(cal[ses_comb == 1], na.rm =T)/sum(cal[],na.rm = T)

sum(cal[ses_comb == 11 | ses_comb == 21], na.rm =T)/1e15
sum(cal[ses_comb == 11 | ses_comb == 21], na.rm =T)/sum(cal[],na.rm = T)

sum(cal[ses_comb == 2], na.rm =T)/1e15
sum(cal[ses_comb == 2], na.rm =T)/sum(cal[],na.rm = T)

sum(cal[ses_comb == 12 | ses_comb == 22], na.rm =T)/1e15
sum(cal[ses_comb == 12 | ses_comb == 22], na.rm =T)/sum(cal[],na.rm = T)

# GDP counts and percentages
sum(gdp[ses_comb == 1], na.rm =T)/1e12
sum(gdp[ses_comb == 1], na.rm =T)/sum(gdp[],na.rm = T)

sum(gdp[ses_comb == 11 | ses_comb == 21], na.rm =T)/1e12
sum(gdp[ses_comb == 11 | ses_comb == 21], na.rm =T)/sum(gdp[],na.rm = T)

sum(gdp[ses_comb == 2], na.rm =T)/1e12
sum(gdp[ses_comb == 2], na.rm =T)/sum(gdp[],na.rm = T)

sum(gdp[ses_comb == 12 | ses_comb == 22], na.rm =T)/1e12
sum(gdp[ses_comb == 12 | ses_comb == 22], na.rm =T)/sum(gdp[],na.rm = T)

# Ramsar counts
sum(ramsar[ses_comb == 1], na.rm =T)
sum(ramsar[ses_comb == 11 | ses_comb == 21], na.rm =T)
sum(ramsar[ses_comb == 2], na.rm =T)
sum(ramsar[ses_comb == 12 | ses_comb == 22], na.rm =T)

# Mean species richness
temp_sr <- amphsr
temp_sr[is.na(ses_comb)] <- NA
temp_sr[ses_comb != 1] <- NA
weighted.mean(x = temp_sr, w = WGS84_areaRaster(0.5), na.rm = T)

temp_sr <- amphsr
temp_sr[is.na(ses_comb)] <- NA
temp_sr[ses_comb != 11 & ses_comb != 21] <- NA
weighted.mean(x = temp_sr, w = WGS84_areaRaster(0.5), na.rm = T)

temp_sr <- amphsr
temp_sr[is.na(ses_comb)] <- NA
temp_sr[ses_comb != 2] <- NA
weighted.mean(x = temp_sr, w = WGS84_areaRaster(0.5), na.rm = T)

temp_sr <- amphsr
temp_sr[is.na(ses_comb)] <- NA
temp_sr[ses_comb != 12 & ses_comb != 22] <- NA
weighted.mean(x = temp_sr, w = WGS84_areaRaster(0.5), na.rm = T)