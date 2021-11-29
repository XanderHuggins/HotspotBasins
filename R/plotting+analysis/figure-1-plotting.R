# Name: figure-1-plotting.R
# Description: Plot three maps: (1) FW stress, (2) Storage trends, (3) Bivariate, and the distributions of social-ecological activity as scatter plots

# Use 'here' package for easy path management
library(here)

# Import all setup and user-defined functions in R/setup and R/udfs folders
invisible(sapply(paste0(here("R/setup"), "/", list.files(here("R/setup"))), source)) 
invisible(sapply(paste0(here("R/udf"), "/", list.files(here("R/udf"))), source))

# Import coastlines and basin shapefile with plotting data
coastlines <- sf::read_sf(here('Data/hybas_l4_coastlines.shp'))
basins_data <- sf::read_sf(here('R/plotting+stats/Basin_data.shp'))

# Create bivariate classes for plotting 
# Reclassify freshwater stress, with class breaks at 10%, 40%
basins_data$fwstrsClass <- rep(NA, nrow(basins_data))
for (i in 1:nrow(basins_data)) { 
  basins_data$fwstrsClass[basins_data$fwstrs < 0.1] <- 0
  basins_data$fwstrsClass[basins_data$fwstrs >= 0.1 & basins_data$fwstrs < 0.4] <- 10
  basins_data$fwstrsClass[basins_data$fwstrs >= 0.4] <- 20
}

# Reclassify storage trend, with class breaks at -3 and 3 mm/yr
basins_data$twsClass <- rep(NA, nrow(basins_data))
for (i in 1:nrow(basins_data)) { 
  basins_data$twsClass[basins_data$tws < -3] <- 0
  basins_data$twsClass[basins_data$tws >= -3 & basins_data$tws < 3] <- 1
  basins_data$twsClass[basins_data$tws >= 3] <- 2
}

basins_data$bivar <- basins_data$fwstrsClass + basins_data$twsClass

# Map of FW stress per basin
tm <- tm_shape(basins_data, projection = "+proj=robin") +
  tm_polygons(col = "fwstrs",  style = "cont", 
              palette = scico(20, palette = "bilbao", direction = 1), 
            midpoint = 0.2, breaks = c(0, 0.4), alpha = 1.00,
            lwd = 0.2, colorNA = 'grey50') +
  tm_shape(coastlines) +
  tm_borders(lwd = 0.7, col = "black") +
  tm_layout(legend.show = F, earth.boundary = c(-179, -60, 179, 88),
            earth.boundary.color = "white", space.color = "white",
            legend.frame = F, frame = F,
            outer.margins = c(-0, -0.09, -0, -0.04)) # B, L, T, R
tm    
tmap_save(tm, here::here("plot_outputs/Fig1_stress.pdf"), dpi = 400, units = "in")

# Map of TWS trend per basin
tm <- tm_shape(basins_data, projection = "+proj=robin") +
  tm_polygons(col = "tws", style = "cont", 
              palette = scico(20, palette = "vikO", direction = -1)[3:18], 
              midpoint = 0, breaks = c(-20, 20), alpha = 1.00, 
              lwd = 0.2) +
  tm_shape(coastlines) +
  tm_borders(lwd = 0.7, col = "black") +
  tm_layout(legend.show = F, earth.boundary = c(-179, -60, 179, 88),
            earth.boundary.color = "white", space.color = "white",
            legend.frame = F, frame = F,
            outer.margins = c(-0, -0.09, -0, -0.04)) # B, L, T, R
tm
tmap_save(tm, here::here("plot_outputs/Fig1_storagetrend.pdf"), dpi = 400, units = "in")

# Map of bivariate combination of TWS and FW stress
# Custom color palette
bivar.pal <- c("#EEDAD9", "#E6E6E6", "#DEE8EB",
               "#BF9290", "#BFBFBF", "#9DB2B9",
               "#7E302E", "#7F7F7F", "#456A76")

tm <- tm_shape(basins_data, projection = "+proj=robin") +
  tm_polygons(col = "bivar", style = "cat", 
              palette = bivar.pal, alpha = 1.00, 
              lwd = 0.2, colorNA = '#FFFFFF') +
  tm_shape(coastlines) +
  tm_borders(lwd = 0.7, col = "black") +
  tm_layout(legend.show = F, earth.boundary = c(-179, -60, 179, 88),
            earth.boundary.color = "white", space.color = "white",
            legend.frame = F, frame = F,
            outer.margins = c(-0, -0.09, -0, -0.04)) # B, L, T, R
tm
tmap_save(tm, here::here("plot_outputs/Fig1_bivar.pdf"), dpi = 400, units = "in")

# Loop through plotting the four dimensions  ---- 
basins_data <- sf::read_sf(here('R/plotting/Basin_data.shp')) %>% 
  as.data.frame() # need in dataframe

pltdims <- c('pop', 'cal', 'gdp', 'amph')
aa = 0.4 # Transparency 
basins_data$Log10Stress <- log10(100*basins_data$fwstrs ) # Convert stress to log10 scale 

for (i in 1:nrow(basins_data)) {
  basins_data$Log10Stress[i] <- max(min(basins_data$Log10Stress[i], 2), 0)
  basins_data$tws[i] <- max(min(basins_data$tws[i], 22), -22)
  
}

basins_data$ramsar[basins_data$ramsar == 0] <- NA # set 0s to NA so they don't appear when plotting
basins_data$amph[basins_data$amph == 0] <- NA # set 0s to NA so they don't appear when plotting

for (i in 1:length(pltdims)) {
  i = 4
  mp <- ggplot(data = basins_data, aes_string(x = 'tws', y = 'Log10Stress', size = pltdims[i] )) +
    annotate("rect", xmin=-Inf, xmax=-3, ymin=1, ymax=2, fill="#7D3232", alpha=aa) +  #upleft
    annotate("rect", xmin=-Inf, xmax=-3, ymin=0, ymax=1, fill="#7D6464", alpha=aa) +  #downleft
    annotate("rect", xmin=3, xmax=Inf, ymin=1,   ymax=2, fill="#19324B", alpha=aa) +  #upright
    annotate("rect", xmin=3, xmax=Inf, ymin=0,   ymax=1, fill="#4B647D", alpha=aa) +  #downright
    annotate("rect", xmin=-3, xmax=3, ymin=1,    ymax=2, fill="#E6E6E6", alpha=aa) +  #centre  
    geom_point(alpha= 0.5, shape = 21, color = 'black', stroke = 2) +
    {if(i != 4)scale_size(range = c(.1, 24))} +
    {if(i == 4)scale_size(range = c(.1, 12))} +
    {if(i == 4)
      geom_point(data = basins_data, aes_string(x = 'tws', y = 'Log10Stress', size = 'ramsar'),
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
         paste0(here::here("plot_outputs/Fig1_"), pltdims[i], ".pdf", sep = ""),
         dpi = 500, width = 6, height = 4, units = "in")
}
