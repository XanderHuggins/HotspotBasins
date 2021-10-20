# FeatureAreaAverage ----
# Calculate area-weighted average value of continuous variable per feature and export as raster 

# Area-weighted feature average ----
FeatureAreaAverage <- function(FT.id, RawDS, AreaDS, operation, varnam) {
  # @FT.id: Input raster with ID values unique to each feature
  # @RawDS: Raw raster dataset to calculate area-weighted averages of
  # @AreaDS: Raster with cell values representing the cell area
  # @operation: Specifies summary function (average or sum)
  # @varnam: Variable name to apply to output raster
  
  temp <- raster::stack(FT.id, RawDS, AreaDS) %>% as.data.frame() %>% 
    set_colnames(c("Feature", "RawData", "Area")) 
  
  if (operation == "mean") {
    Sumtab <- temp %>%
      group_by(Feature) %>%
      summarise(
        Stat1 = weighted.mean(x = RawData, w = Area, na.rm = T)
      )  
  }
  
  if (operation == "sum") {
    Sumtab <- temp %>%
      group_by(Feature) %>%
      summarise(
        Stat1 = sum(x = RawData, w = Area, na.rm = T)
      )  
  }
  
  # Format for raster reclassification
  Sumtab$uplim <- Sumtab$Feature+0.5
  Sumtab$lowlim <- Sumtab$Feature-0.5
  rclmtx <- matrix(c(Sumtab$lowlim, Sumtab$uplim, Sumtab$Stat1), ncol = 3)
  
  # Reclassify with summary statistic
  ProdRas <- reclassify(FT.id, rclmtx)
  ProdRas[is.na(FT.id)] <- NA
  names(ProdRas) <- as.character(varnam)
  plot(ProdRas)
  
  # Return raster
  return(ProdRas)
  
}