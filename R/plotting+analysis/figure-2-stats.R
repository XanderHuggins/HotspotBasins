# Name: figure-2-stats.R
# Description: Calculate summary statistics associated with Figure 2

# Use 'here' package for easy path management
library(here)

# Import all setup and user-defined functions in R/setup and R/udfs folders
invisible(sapply(paste0(here("R/setup"), "/", list.files(here("R/setup"))), source)) 
invisible(sapply(paste0(here("R/udf"), "/", list.files(here("R/udf"))), source))

# Import basin shapefile with plotting data
basins_data <- sf::read_sf(here('R/plotting+stats/Basin_data.shp')) %>% 
  as.data.frame()
basins_data <- basins_data[complete.cases(basins_data$fw_stts),]

# Reproduce bivairate classes
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

table(basins_data$bivar) # Provides summary of number of basins in each bivariate class

basins_data$corner <- rep(NA)
basins_data$corner[basins_data$bivar == 1] <- 'LM'
basins_data$corner[basins_data$bivar == 11 | basins_data$bivar == 21] <- 'UM'
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
