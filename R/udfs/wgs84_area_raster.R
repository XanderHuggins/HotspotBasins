# WGS84_areaRaster ----
# Generate raster with cell values representing grid areas based on WGS84 reference ellipsoid. 
# This method is slightly more accurate & sophisticated than that provided by raster::area() function.
# Follows methods outlined in Santini et al. https://doi.org/10.1111/j.1467-9671.2010.01200.x

WGS84_areaRaster <- function(ResT) {
  # @ResT: Desired resolution
  
  library(tidyverse)
  require(raster)
  
  pi.g <- 3.14159265358979  
  Fl <- 0.00335281066474 # Flattening
  SMA <- 6378137.0 # Semi major axis
  e <- sqrt((2*Fl) - (Fl^2)) # Eccentricity  
  RES <- ResT
  
  # error check entry
  if ((90/RES) %% 1 > 1e-6) { stop("'ResT' must a factor of 90") }
  # if (90/RES > (90/(1/24))) { stop("'ResT' is too fine (will require more memory than can be allocated") }
  
  # initialize dataframe with geodetic latitudes
  df_a <- data.frame(LowLAT = seq(-90, 90-RES, by = RES), 
                     UppLAT = seq(-90+RES, 90, by = RES))
  
  # Convert geodetic latitudes degrees to radians
  df_a$LowLATrad <- df_a$LowLAT * pi.g / 180
  df_a$UppLATrad <- df_a$UppLAT * pi.g / 180
  
  # Calculate q1 and q2
  df_a$q1 <- (1-e*e)* ((sin(df_a$LowLATrad)/(1-e*e*sin(df_a$LowLATrad)^2)) - ((1/(2*e))*log((1-e*sin(df_a$LowLATrad))/(1+e*sin(df_a$LowLATrad)))))
  
  df_a$q2 <- (1-e*e)* ((sin(df_a$UppLATrad)/(1-e*e*sin(df_a$UppLATrad)^2)) - ((1/(2*e))*log((1-e*sin(df_a$UppLATrad))/(1+e*sin(df_a$UppLATrad)))))
  
  # calculate q constant
  q_const <- (1-e*e)* ((sin(pi.g/2)/(1-e*e*sin(pi.g/2)^2)) - ((1/(2*e))*log((1-e*sin(pi.g/2))/(1 + e*sin(pi.g/2)))))
  
  # Calculate authaltic latitudes
  df_a$phi1 <- asin(df_a$q1 / q_const)
  df_a$phi2 <- asin(df_a$q2 / q_const)
  
  # Calculate authaltic radius
  R.adius <- sqrt(SMA*SMA*q_const/2)
  
  # Calculate cell size in radians
  CS <- (RES) * pi.g/180
  
  # Calculate cell area in m2
  df_a$area_m2 <- R.adius*R.adius*CS*(sin(df_a$phi2)-sin(df_a$phi1))
  
  # Convert to raster, and replicate column throughout global longitude domain
  WGS84area_km2 <- matrix(df_a$area_m2/1e6, nrow = 180/RES, ncol = 360/RES, 
                          byrow = FALSE, dimnames = NULL) %>% raster()
  extent(WGS84area_km2) <- c(-180, 180, -90, 90) # set extent of raster
  crs(WGS84area_km2) <- crs("+proj=longlat") # set CRS of raster
  
  WGS84area_km2 <- raster::flip(WGS84area_km2, direction = "y")
  
  message(paste0("Calculated global surface area at: ", RES, 
                 "deg. is ", sum(WGS84area_km2[]), " km2.", sep = ""))
  
  return(WGS84area_km2)
}