# Name: figure-3-plotting.R
# Description: Plot three maps: (1) Hotspot basins map, (2) Ecological vulnerability, (3) Social vulnerability

# Use 'here' package for easy path management  
library(here)

# Import all setup and user-defined functions in R/setup and R/udfs folders
invisible(sapply(paste0(here("R/setup"), "/", list.files(here("R/setup"))), source)) 
invisible(sapply(paste0(here("R/udf"), "/", list.files(here("R/udf"))), source))

# Import coastlines and basin shapefile with plotting data
coastlines <- sf::read_sf(here('Data/hybas_l4_coastlines.shp'))
basins_data <- sf::read_sf(here('R/plotting+stats/Basin_data.shp'))

# Custom palette
cus.pal <- c("grey80", "#B99232", "#9D5819", "#7E1900")

# Plot hotspot basins map
basins_data$hotspots <- rep(NA)
basins_data$hotspots[basins_data$ss_vcls >= 2] <- 1
basins_data_m <- basins_data %>% filter(hotspots == 1)

# dissolve for outlining
basins_data_m <- fasterize(sf = basins_data_m, raster = WGS84_areaRaster(0.5))
basins_data_m <- raster::rasterToPolygons(basins_data_m, na.rm = T, dissolve = T)

tm <- tm_shape(basins_data, projection = "+proj=robin") +
  tm_polygons(col = "ss_vcls", style = "cat", palette = cus.pal,
              lwd = 0.2, colorNA = '#FFFFFF') +
  tm_shape(basins_data_m) + tm_borders(lwd = 2, col = 'black') + 
  tm_shape(coastlines) +
  tm_borders(lwd = 0.7, col = "black") +
  tm_layout(legend.show = F, earth.boundary = c(-179, -60, 179, 88),
            earth.boundary.color = "white", space.color = "white",
            legend.frame = F, frame = F,
            outer.margins = c(-0, -0.09, -0, -0.04)) # B, L, T, R
tm
tmap_save(tm, here::here("plot_outputs/Fig3_hotspotbasins.pdf"), dpi = 400, units = "in")

# Plot ecological vulnerability classes
tm <- tm_shape(basins_data, projection = "+proj=robin") +
  tm_polygons(col = "ecl_vcl", style = "cat", palette = cus.pal,
              lwd = 0.2, colorNA = '#FFFFFF') +
  tm_shape(coastlines) +
  tm_borders(lwd = 0.7, col = "black") +
  tm_layout(legend.show = F, earth.boundary = c(-179, -60, 179, 88),
            earth.boundary.color = "white", space.color = "white",
            legend.frame = F, frame = F,
            outer.margins = c(-0, -0.09, -0, -0.04)) # B, L, T, R
tm
tmap_save(tm, here::here("plot_outputs/Fig3_ecol_vuln.pdf"), dpi = 400, units = "in")

# Plot social vulnerability classes
tm <- tm_shape(basins_data, projection = "+proj=robin") +
  tm_polygons(col = "scl_vcl", style = "cat", palette = cus.pal,
              lwd = 0.2, colorNA = '#FFFFFF') +
  tm_shape(coastlines) +
  tm_borders(lwd = 0.7, col = "black") +
  tm_layout(legend.show = F, earth.boundary = c(-179, -60, 179, 88),
            earth.boundary.color = "white", space.color = "white",
            legend.frame = F, frame = F,
            outer.margins = c(-0, -0.09, -0, -0.04)) # B, L, T, R
tm
tmap_save(tm, here::here("plot_outputs/Fig3_socl_vuln.pdf"), dpi = 400, units = "in")


# Make ribbon plot on classification ----
basins_data <- basins_data[complete.cases(basins_data$fw_stts), ]

# Need to re-identify head/tail break values
htb_ses <- ht_breaks(x = basins_data$ses_vln, tsh = 0.8)

# Generate points along curves for vulnerability class thresholds
invlin <- data.frame(x = seq(htb_ses[1], 1, length.out = 1000), 
                     y1 = rep(NA, 1000),
                     y2 = rep(NA, 1000),
                     y3 = rep(NA, 1000))
for (i in 1:1000) { 
  invlin$y1[i] <- min(htb_ses[1]/invlin$x[i], 1)
  invlin$y2[i] <- min(htb_ses[2]/invlin$x[i], 1)
  invlin$y3[i] <- min(htb_ses[3]/invlin$x[i], 1)
}

tsh_ln <- invlin[invlin[,1] >= htb_ses[2]-0.0025,]

aa = 0.4
mp <- ggplot() +
  geom_ribbon(data = invlin, aes(x = x, ymin = y1, ymax = y2), fill = "#B99232", alpha = aa) +
  geom_ribbon(data = invlin, aes(x = x, ymin = y2, ymax = y3), fill = "#9D5819", alpha = aa) +
  geom_ribbon(data = invlin, aes(x = x, ymin = y3, ymax = 1),  fill = "#7E1900", alpha = aa) +
  geom_point(data = basins_data, aes(x = fw_stts, y = ses_sns, fill = as.factor(ss_vcls), 
                                     color = as.factor(ss_vcls)), 
             alpha = 1, shape = 21, size = 3, stroke = 0.5) +
  scale_fill_manual(values = c('#CCCCCC', '#B99232', '#9D5819', '#7E1900')) +
  scale_color_manual(values = c('#666666', '#666666', 'black', 'black')) +
  geom_point(data = subset(basins_data, ss_vcls >= 2), aes(x = fw_stts, y = ses_sns),
             alpha= 1, shape = 21, color = 'black', size = 3, stroke = 1.5) +
  geom_line(data = invlin, aes(x = x, y = y1), col = 'black', lwd = 0.5, linetype = "dashed") +
  geom_line(data = tsh_ln, aes(x = x, y = y2), col = 'black', lwd = 1.5) +
  geom_line(data = invlin, aes(x = x, y = y3), col = 'black', lwd = 0.5, linetype = "dashed") +
  scale_size(range = c(.1, 24)) +
  coord_cartesian(xlim = c(0, 1.05), ylim = c(0, 1.05), expand = c(0, 0), clip = "off") +
  scale_y_continuous(breaks = seq(0, 1, 0.2), 
                     limits = c(0, 1)) + 
  theme1 + theme(axis.ticks.x = element_line(size = 1)) 
mp
ggsave(plot = last_plot(), 
       here::here("plot_outputs/Fig3_ribbon.pdf"),
       dpi = 500, width = 6, height = 4, units = "in")


# Make histograms for basin counts and surface area
s_df <- basins_data %>% as.data.frame() %>% 
  dplyr::group_by(ss_vcls) %>% 
  summarize(
    b_count = unique(pfaf_id, na.rm = T) %>% length(),
    s_area = sum(area, na.rm = T)
  )
s_df$sa_pct <- s_df$s_area/sum(s_df$s_area)
head(s_df)

ggplot(data = s_df, aes(x = as.factor(ss_vcls), y = b_count, fill = as.factor(ss_vcls))) +
  geom_bar(stat = 'identity', color = 'black', lwd = 3.5, width = 1) +
  scale_fill_manual(values = c('#CCCCCC', '#B99232', '#9D5819', '#7E1900')) +
  coord_cartesian(expand = c(0, 0), clip = "off") +
  theme1 + theme(axis.ticks.x = element_line(size = 1)) 
ggsave(plot = last_plot(), 
       here::here("plot_outputs/Fig3_hist_basincount.pdf"),
       dpi = 500, width = 6, height = 4, units = "in")

ggplot(data = s_df, aes(x = as.factor(ss_vcls), y = s_area, fill = as.factor(ss_vcls))) +
  geom_bar(stat = 'identity', color = 'black', lwd = 3.5, width = 1) +
  scale_fill_manual(values = c('#CCCCCC', '#B99232', '#9D5819', '#7E1900')) +
  coord_cartesian(expand = c(0, 0), clip = "off") +
  theme1 + theme(axis.ticks.x = element_line(size = 1)) 
ggsave(plot = last_plot(), 
       here::here("plot_outputs/Fig3_hist_surfacearea.pdf"),
       dpi = 500, width = 6, height = 4, units = "in")
