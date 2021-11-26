# Name: figure-4-stats.R
# Description: Calculate summary statistics associated with Figure 4

# Use 'here' package for easy path management
library(here)

# Import all setup and user-defined functions in R/setup and R/udfs folders
invisible(sapply(paste0(here("R/setup"), "/", list.files(here("R/setup"))), source)) 
invisible(sapply(paste0(here("R/udf"), "/", list.files(here("R/udf"))), source))

# Import basin shapefile with plotting data
basins_data <- sf::read_sf(here('R/plotting+stats/Basin_iwrm_data.shp')) %>% 
  as.data.frame()
basins_data <- basins_data[complete.cases(basins_data$ss_vcls),]

htb_ses <- ht_breaks(x = basins_data$ses_vln, tsh = 0.8)

# Identify number of basins 
basins_data %>% nrow() ## output = 1204; total number of basins
basins_data %>% filter(ss_vcls >= 2) %>% nrow() ## output = 168; number of hotspot basins
basins_data %>% filter(ss_vcls >= 2 & tsbd >= 2) %>% nrow() ## output = 61; # number of hotspot basins that are transboundary
61/168
# output = 36%

# Calculate average IWRM score for non-transboundary hotspot basins
basins_data %>% filter(ss_vcls >= 2 & tsbd == 1) %>% pull(iwrm) %>% mean() # output = 55.5
basins_data %>% filter(ss_vcls >= 2 & tsbd > 1) %>% pull(iwrm) %>% mean() # output = 50.3 
