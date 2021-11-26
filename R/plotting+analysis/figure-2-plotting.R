# Name: figure-2-plotting.R
# Description: Plot three maps: (1) Adaptive capacity, (2) Basin freshwater status, (3) Bivariate, and the distributions of social-ecological activity as scatter plots

# Use 'here' package for easy path management  
library(here)

# Import all setup and user-defined functions in R/setup and R/udfs folders
invisible(sapply(paste0(here("R/setup"), "/", list.files(here("R/setup"))), source)) 
invisible(sapply(paste0(here("R/udf"), "/", list.files(here("R/udf"))), source))

# Import coastlines and basin shapefile with plotting data
coastlines <- sf::read_sf(here('Data/hybas_l4_coastlines.shp'))
basins_data <- sf::read_sf(here('R/plotting+stats/Basin_data.shp'))

# To create bivariate legend, we use percentiles as break values
# Create dataframe to hold these break values
# Calculate percentile values
ptls <- data.frame(P = c(33, 67, 80), 
                   ac  = rep(NA, 3),
                   status = rep(NA, 3))

for (i in 1:nrow(ptls)) {
  ptls$ac[i]  <-    quantile(x = basins_data$ac,   probs = ptls$P[i]/100, na.rm = T)
  ptls$status[i] <- quantile(x = basins_data$fw_stts, probs = ptls$P[i]/100, na.rm = T)
}

# Reclassify adaptive capacity with class breaks at P33 and P67 
basins_data$acClass <- rep(NA, nrow(basins_data))
for (i in 1:nrow(basins_data)) { 
  basins_data$acClass[basins_data$ac < ptls$ac[1]] <- 0
  basins_data$acClass[basins_data$ac >= ptls$ac[1] & basins_data$ac < ptls$ac[2]] <- 10
  basins_data$acClass[basins_data$ac >= ptls$ac[2]] <- 20
}

# Reclassify basin freshwater status with class breaks at P67 and P80
basins_data$statusClass <- rep(NA, nrow(basins_data))
for (i in 1:nrow(basins_data)) { 
  basins_data$statusClass[basins_data$fw_stts < ptls$status[2]] <- 0
  basins_data$statusClass[basins_data$fw_stts >= ptls$status[2] & basins_data$fw_stts < ptls$status[3]] <- 1
  basins_data$statusClass[basins_data$fw_stts >= ptls$status[3]] <- 2
}

basins_data$bivar <- basins_data$acClass + basins_data$statusClass

# Map of inverted adaptive capacity
basins_data$inv_ac <- 1 - basins_data$ac

tm <- tm_shape(basins_data, projection = "+proj=robin") +
  tm_polygons(col = "inv_ac",  style = "cont", 
              palette = scico(20, palette = "batlow", direction = 1),
              breaks= c(0,1), midpoint = 0.5,
              lwd = 0.2, colorNA = '#b3b3b3') +
  tm_shape(coastlines) +
  tm_borders(lwd = 0.7, col = "black") +
  tm_layout(legend.show = F, earth.boundary = c(-179, -60, 179, 88),
            earth.boundary.color = "white", space.color = "white",
            legend.frame = F, frame = F,
            outer.margins = c(-0, -0.09, -0, -0.04)) # B, L, T, R
tm
tmap_save(tm, here::here("plot_outputs/Fig2_ac.pdf"), dpi = 400, units = "in")


# Basin freshwater status
basins_data$mask <- rep(NA)
basins_data$mask[basins_data$fw_stts < 0.05] <- 1
basins_data_m <- basins_data %>% filter(mask == 1)

tm <- tm_shape(basins_data, projection = "+proj=robin") +
  tm_polygons(col = "fw_stts",  style = "cont", 
              palette = scico(20, palette = "lajolla", direction = 1)[1:17], 
              midpoint = 0.5, breaks = c(0, 1), alpha = 1.00,
              lwd = 0.2) +
  tm_shape(basins_data_m) + tm_polygons(col = '#b3b3b3', lwd = 0.2) + 
  tm_shape(coastlines) +
  tm_borders(lwd = 0.7, col = "black") +
  tm_layout(legend.show = F, earth.boundary = c(-179, -60, 179, 88),
            earth.boundary.color = "white", space.color = "white",
            legend.frame = F, frame = F,
            outer.margins = c(-0, -0.09, -0, -0.04)) # B, L, T, R
tm
tmap_save(tm, here::here("plot_outputs/Fig2_fw_status.pdf"), dpi = 400, units = "in")


# Bivariate plot
bivar.pal <- c("#7F7F7F", "#19255A", "#7E1900",# bottom row of bivar legend, L to R
             "#BFBFBF", "#455186", "#9D5819",
             "#E6E6E6", "#7683B5", "#B99232")# upper row of bivar legend, L to R

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
tmap_save(tm, here::here("plot_outputs/Fig2_bivar.pdf"), dpi = 400, units = "in")


# Loop through plotting the four dimensions  ---- 
basins_data <- basins_data %>% 
  as.data.frame() # need in dataframe
basins_data <- basins_data[complete.cases(basins_data$fw_stts),]

pltdims <- c('pop', 'cal', 'gdp', 'amph')
aa = 0.4 # Transparency 

basins_data$ramsar[basins_data$ramsar == 0] <- NA # set 0s to NA so they don't appear when plotting
basins_data$amph[basins_data$amph == 0] <- NA # set 0s to NA so they don't appear when plotting

for (i in 1:length(pltdims)) {
  
  mp <- ggplot(data = basins_data, aes_string(x = 'fw_stts', y = 'ac', size = pltdims[i] )) +
    annotate("rect", xmin=0,            xmax=ptls$status[2],   ymin=ptls$ac[2], ymax=1,           fill="#E6E6E6", alpha=aa) +  
    annotate("rect", xmin=0,            xmax=ptls$status[2],   ymin=ptls$ac[1], ymax=ptls$ac[2], fill="#BFBFBF", alpha=aa) +  
    annotate("rect", xmin=0,            xmax=ptls$status[2],   ymin=0,           ymax=ptls$ac[1], fill="#7F7F7F", alpha=aa) +     
    annotate("rect", xmin=ptls$status[2], xmax=ptls$status[3], ymin=ptls$ac[2], ymax=1,           fill="#7683B5", alpha=aa) +     
    annotate("rect", xmin=ptls$status[2], xmax=ptls$status[3], ymin=ptls$ac[1], ymax=ptls$ac[2], fill="#455186", alpha=aa) +
    annotate("rect", xmin=ptls$status[2], xmax=ptls$status[3], ymin=0,           ymax=ptls$ac[1], fill="#19255A", alpha=aa) +
    annotate("rect", xmin=ptls$status[3], xmax=1,              ymin=ptls$ac[2], ymax=1,           fill="#B99232", alpha=aa) +  
    annotate("rect", xmin=ptls$status[3], xmax=1,              ymin=ptls$ac[1], ymax=ptls$ac[2], fill="#9D5819", alpha=aa) +  
    annotate("rect", xmin=ptls$status[3], xmax=1,              ymin=0,           ymax=ptls$ac[1], fill="#7E1900", alpha=aa) +
    geom_point(alpha= 0.5, shape = 21, color = 'black', stroke = 2) +
    {if(i != 4)scale_size(range = c(0.1, 24))} +
    {if(i == 4)scale_size(range = c(0.1, 12))} +
    {if(i == 4)
      geom_point(data = basins_data, aes_string(x = 'fw_stts', y = 'ac', size = 'ramsar'),
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
         paste0(here::here("plot_outputs/Fig2_"), pltdims[i], ".pdf", sep = ""),
         dpi = 500, width = 6, height = 4, units = "in")
}
