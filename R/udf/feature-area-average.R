# Name: FeatureAreaAverage 
# Description: Calculates the area-weighted average or sum of a continuous raster based on another raster specifying zonal unit IDs, and exports a raster with the specified summary statistic

FeatureAreaAverage <- function(FT.id, RawDS, AreaDS, operation, varnam) {
  # Function arguments:
  # FT.id: Raster with values specifying zonal unit unique IDs
  # RawDS: Raster to summarize
  # AreaDS: Raster with cell values representing cell-specific surface area
  # operation: Choice between "average" (which is always area-weighted) or  "sum"
  # varnam: Variable name to set for the output summary statistic raster raster
  
  # Stack the feature ID, raster to summarize, and area raster, and convert to dataframe
  temp <- raster::stack(FT.id, RawDS, AreaDS) %>% as.data.frame() %>% 
    set_colnames(c("Feature", "RawData", "Area")) 
  
  # Calculate area-weighted mean per feature if specified
  if (operation == "mean") {
    Sumtab <- temp %>%
      group_by(Feature) %>%
      summarise(
        Stat1 = weighted.mean(x = RawData, w = Area, na.rm = T)
      )  
  }
  
  # Calculate sum per feature if specified
  if (operation == "sum") {
    Sumtab <- temp %>%
      group_by(Feature) %>%
      summarise(
        Stat1 = sum(x = RawData, na.rm = T)
      )  
  }
  
  # Create a reclassification matrix to convert dataframe results to raster
  Sumtab$uplim <- Sumtab$Feature+0.5
  Sumtab$lowlim <- Sumtab$Feature-0.5
  rclmtx <- matrix(c(Sumtab$lowlim, Sumtab$uplim, Sumtab$Stat1), ncol = 3)
  
  # Reclassify feature ID raster with summary statistic results
  ProdRas <- reclassify(FT.id, rclmtx)
  ProdRas[is.na(FT.id)] <- NA
  names(ProdRas) <- as.character(varnam)
  plot(ProdRas)
  
  # Return raster
  return(ProdRas)
  
}
