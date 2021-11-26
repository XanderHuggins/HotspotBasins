# Name: figure-4-plotting.R
# Description: Plot a global map of IWRM implementation and hotspot basins, and a scatterplot comparing the two at a basin scale

# Use 'here' package for easy path management
library(here)

# Import all setup and user-defined functions in R/setup and R/udfs folders
invisible(sapply(paste0(here("R/setup"), "/", list.files(here("R/setup"))), source)) 
invisible(sapply(paste0(here("R/udf"), "/", list.files(here("R/udf"))), source))

# Import coastlines and basin shapefile with plotting data
coastlines <- sf::read_sf(here('Data/hybas_l4_coastlines.shp'))
basins_data <- sf::read_sf(here('R/plotting+stats/Basin_iwrm_data.shp'))
iwrm_harm <- raster(here::here('Data/iwrm_fill.tif'))

# Classify IWRM based on IWRM data portal classes
rclmat <- c(0,  10, 1,  
            10, 30, 2,
            30, 50, 3,
            50, 70, 4,
            70, 90, 5,
            90, 100,6) %>% matrix(ncol = 3, byrow = T)
IWRMclass <- reclassify(iwrm_harm, rclmat)

IWRMclass <- rasterToPolygons(IWRMclass, dissolve = T)

# Identify non-hotspots and shade on map
basins_data$faderegs <- rep(NA)
basins_data$faderegs[basins_data$ss_vcls < 2] <- 1
basins_data_d <- basins_data %>% filter(faderegs == 1)

# Identify non-hotspots and shade on map
basins_data$hotspots <- rep(NA)
basins_data$hotspots[basins_data$ss_vcls >= 2] <- 1
basins_data_h <- basins_data %>% filter(hotspots == 1)
basins_data_h <- fasterize(sf = basins_data_h, raster = WGS84_areaRaster(0.5))
basins_data_h <- raster::rasterToPolygons(basins_data_h, na.rm = T, dissolve = T)

# Custom palette
cus.pal <- c("#D57F31", "#E8BA24", "#81BF81", "#6CBDEF", "#0D7EB6")

tm <- tm_shape(coastlines, projection = "+proj=robin") + tm_fill(col = '#b3b3b3') +
  tm_shape(IWRMclass) +
  tm_polygons(col = "iwrm_fill", style = "cat", palette = cus.pal,
              lwd = 0, colorNA = '#b3b3b3') +
  tm_shape(basins_data) + tm_borders(lwd = 0.2) + 
  tm_shape(basins_data_d) + tm_fill(col = 'black', alpha = 0.5) + 
  tm_shape(basins_data_h) + tm_borders(col = 'black', lwd = 2) + 
  tm_shape(coastlines) +
  tm_borders(lwd = 0.7, col = "black") +
  tm_layout(legend.show = F, earth.boundary = c(-179, -60, 179, 88),
            earth.boundary.color = "white", space.color = "white",
            legend.frame = F, frame = F,
            outer.margins = c(-0, -0.09, -0, -0.04)) # B, L, T, R
tm
tmap_save(tm, here::here("plot_outputs/Fig4_iwrm-map.pdf"), dpi = 400, units = "in")


# Make scatter plot ---- 

# Need class breaks 
basins_data_c <- basins_data[complete.cases(basins_data$ss_vcls),]
htb_ses <- ht_breaks(x = basins_data_c$ses_vln, tsh = 0.8)

aa = 0.4

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
  geom_vline(xintercept = htb_ses[1], lwd = 1, linetype = 'dashed') +
  geom_vline(xintercept = htb_ses[2], lwd = 1, linetype = 'dashed') +
  geom_vline(xintercept = htb_ses[3], lwd = 1, linetype = 'dashed') +
  # add red outline to transboundary basins
  geom_point(data = subset(basins_data, tsbd >= 2), aes(x = ses_vln, y = iwrm), 
             alpha= 1, shape = 21, color = 'firebrick3', size = 3, lwd = 4) +
  
  geom_point(data = subset(basins_data, tsbd >= 3), aes(x = ses_vln, y = iwrm), 
             alpha= 1, shape = 21, color = 'firebrick3', size = 4, lwd = 4) +
  
  geom_point(data = subset(basins_data, tsbd >= 4), aes(x = ses_vln, y = iwrm), 
             alpha= 1, shape = 21, color = 'firebrick3', size = 5, lwd = 4) +
  
  geom_point(data = subset(basins_data, tsbd >= 5), aes(x = ses_vln, y = iwrm), 
             alpha= 1, shape = 21, color = 'firebrick3', size = 6, lwd = 4) +
  
  geom_point(data = subset(basins_data, tsbd >= 6), aes(x = ses_vln, y = iwrm), 
             alpha= 1, shape = 21, color = 'firebrick3', size = 7, lwd = 4) +
  # plot all basins satisfying mask constraint 
  
  # semi-opaque for non-transboundary
  geom_point(data = subset(basins_data, tsbd == 1), aes(x = ses_vln, y = iwrm), 
             alpha= 0.25, shape = 19, color = 'black', size = 2) +
  
  # full opacity for transboundary
  geom_point(data = subset(basins_data, tsbd >= 2), aes(x = ses_vln, y = iwrm), 
             alpha= 1, shape = 19, color = 'black', size = 2) +
  
  # mask areas under first break
  annotate("rect", xmin=0,  xmax=htb_ses[2], ymin=0,  ymax=100, fill="white", alpha=1.5*aa) +
  scale_y_continuous(breaks = seq(10, 90, 20), 
                     limits = c(0, 100)) + 
  theme1 + theme(axis.ticks.x = element_line(size = 1)) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 100), expand = c(0, 0), clip = "on")
mp
ggsave(plot = mp, 
       here::here("plot_outputs/Fig4_iwrm-scatter.pdf"),
       bg = 'transparent', dpi = 500, width = 6, height = 4, units = "in")
