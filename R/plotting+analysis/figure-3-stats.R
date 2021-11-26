# Name: figure-3-stats.R
# Description: Calculate summary statistics associated with Figure 3

# Use 'here' package for easy path management
library(here)

# Import all setup and user-defined functions in R/setup and R/udfs folders
invisible(sapply(paste0(here("R/setup"), "/", list.files(here("R/setup"))), source)) 
invisible(sapply(paste0(here("R/udf"), "/", list.files(here("R/udf"))), source))

# Import basin shapefile with plotting data
basins_data <- sf::read_sf(here('R/plotting+stats/Basin_data.shp')) %>% 
  as.data.frame()

sum_stats <- basins_data %>% 
  group_by(ss_vcls) %>% 
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

# Extract Ramsar Sites within hotspot basins
basins_data <- sf::read_sf(here('R/plotting+stats/Basin_data.shp'))
hotspots <- fasterize(basins_data, raster = WGS84_areaRaster(0.5), field = 'ss_vcls')

ramsar <- readr::read_csv('D:/!! Geodatabase/Biodiversity-Ecosystems/Ramsar-wetlands/ris-tabular-20210322.csv')
ramsar.pts <- SpatialPoints(cbind(ramsar$Longitude, ramsar$Latitude))
crs(ramsar.pts) <- crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+")

ExtValue <- extract(hotspots, ramsar.pts)

RamsarHotspots <- cbind(as.data.frame(ramsar),ExtValue) 
RamsarHotspots <- RamsarHotspots %>% 
  as.data.frame() %>% 
  filter(ExtValue >= 2) %>% 
  dplyr::select('Ramsar Site No.', 'Site name', 'Latitude', 'Longitude', 'Wetland Type', 'ExtValue')
names(RamsarHotspots)[6] <- 'Vulnerability_class_n'

RamsarHotspots$Vulnerability_class_t <- rep(NA)
RamsarHotspots$Vulnerability_class_t[RamsarHotspots$Vulnerability_class_n == 3] <- "Hotspot; Very high vulnerability"
RamsarHotspots$Vulnerability_class_t[RamsarHotspots$Vulnerability_class_n == 2] <- "Hotspot; High vulnerability"

write.csv(RamsarHotspots, here('plot_outputs/RamsarHotspots.csv'))
