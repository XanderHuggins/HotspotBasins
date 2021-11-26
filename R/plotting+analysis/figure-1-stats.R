# Name: figure-1-stats.R
# Description: Calculate summary statistics associated with Figure 1

# Use 'here' package for easy path management
library(here)

# Import all setup and user-defined functions in R/setup and R/udfs folders
invisible(sapply(paste0(here("R/setup"), "/", list.files(here("R/setup"))), source)) 
invisible(sapply(paste0(here("R/udf"), "/", list.files(here("R/udf"))), source))

# Import basin shapefile with plotting data
basins_data <- sf::read_sf(here('R/plotting+stats/Basin_data.shp')) %>% 
  as.data.frame()

# Reproduce bivairate classes
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

table(basins_data$bivar) # Provides summary of number of basins in each bivariate class

basins_data$corner <- rep(NA)
basins_data$corner[basins_data$bivar == 0] <- 'LL'
basins_data$corner[basins_data$bivar == 10 | basins_data$bivar == 20] <- 'UL'
basins_data$corner[basins_data$bivar == 2] <- 'LR'
basins_data$corner[basins_data$bivar == 12 | basins_data$bivar == 22] <- 'UR'

sum_stats <- basins_data %>% 
  group_by(corner) %>% 
  summarise(
    pop = sum(pop, na.rm = T),
    popPCT = sum(pop, na.rm = T)/sum(basins_data$pop, na.rm = T),
    cal = sum(cal, na.rm = T),
    calPCT = sum(cal, na.rm = T)/sum(basins_data$cal, na.rm = T),
    gdp = sum(gdp, na.rm = T),
    gdpPCT = sum(gdp, na.rm = T)/sum(basins_data$gdp, na.rm = T),
    amph = weighted.mean(amph, w = area, na.rm = T),
    rmsr = sum(ramsar, na.rm = T)
  )

sum_stats
