# RasterAreaPercentiles ----
# Reclassify a raster into area-weighted percentiles 

RasterAreaPercentiles <- function(RasterToClassify, WeightRaster, MaskRaster, clipToExtent, CRS.set, ext.set){
  # @RasterToClassify: Input raster
  # @WeightRaster: Typically representing cell area
  # @MaskRaster: Binary raster indicating cells to mask
  # @clipToExtent: If == "clip", use MaskRaster
  # @CRS.set: Coordinate reference system
  # @ext.set: Optional, if wanting to specify the extent of the output raster
  
  PercentileRaster <- raster(RasterToClassify) # Initialize percentile raster
  
  crs(RasterToClassify) <- CRS.set
  extent(RasterToClassify) <- ext.set
  message(paste0("CRS and extent set to: ", crs(RasterToClassify), "  &  ",
                 extent(RasterToClassify), sep = ""))
  
  RasterToClassify[MaskRaster != 1] <- NA
  
  m.df <- raster::stack(RasterToClassify, WeightRaster, MaskRaster) %>% 
    as.data.frame() %>% 
    set_colnames(c("Input", "Weight", "Mask"))
  m.df <- m.df[complete.cases(m.df$Input),]
  
  if(clipToExtent == "clip"){
    m.df <- m.df %>% dplyr::filter(Mask == 1)
  }
  
  pb <- txtProgressBar(min = 0, max = 99, style = 3)
  
  for(i in 0:99){
    j = 1 - (i*0.01)
    k = 0.99 - (i*0.01)
    
    # Upper bound
    ub = as.numeric(unname(spatstat::weighted.quantile(m.df$Input, 
                                                       m.df$Weight, 
                                                       j, 
                                                       na.rm = TRUE)))
    
    # Lower bound
    lb = as.numeric(unname(spatstat::weighted.quantile(m.df$Input, 
                                                       m.df$Weight, 
                                                       k, 
                                                       na.rm = TRUE)))
    
    PercentileRaster[RasterToClassify <= ub & RasterToClassify > lb] <- j
    setTxtProgressBar(pb, i)
  }
  
  PercentileRaster[is.na(PercentileRaster)] <- 0
  PercentileRaster[is.na(MaskRaster) | MaskRaster != 1] <- NA # mask classified by mask raster
  plot(PercentileRaster)
  
  return(PercentileRaster)
}