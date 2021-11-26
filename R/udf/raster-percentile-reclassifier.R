# Name: RasterAreaPercentiles 
# Description: Reclassifies a raster based on its area-weighted percentiles 

RasterAreaPercentiles <- function(RasterToClassify, WeightRaster, MaskRaster, clipToExtent = NULL, CRS.set, ext.set){
  # Function arguments:
  # RasterToClassify: Input raster to reclassify
  # WeightRaster: Raster with values representing cell area
  # MaskRaster: Binary raster representing cells to mask from function output (values == 1 are retained)
  # clipToExtent: If == "clip", use MaskRaster
  # CRS.set: Coordinate reference system
  # ext.set: Optional, if wanting to specify the extent of the output raster
  
  PercentileRaster <- raster(RasterToClassify) # Initialize percentile raster
  crs(RasterToClassify) <- CRS.set # Set CRS
  extent(RasterToClassify) <- ext.set # Set extent
  message(paste0("CRS and extent set to: ", crs(RasterToClassify), "  &  ",
                 extent(RasterToClassify), sep = ""))
  
  RasterToClassify[MaskRaster != 1] <- NA # Apply mask
  
  # Stack rasters and convert to dataframe for faster processing
  m.df <- raster::stack(RasterToClassify, WeightRaster, MaskRaster) %>% 
    as.data.frame() %>% 
    set_colnames(c("Input", "Weight", "Mask"))
  m.df <- m.df[complete.cases(m.df$Input),]
  
  # Remove cells outside of mask domain from percentile calculations
  if(clipToExtent == "clip"){
    m.df <- m.df %>% dplyr::filter(Mask == 1)
  }
  
  pb <- txtProgressBar(min = 0, max = 99, style = 3)
  
  # Loop through each percentile
  for(i in 0:99){
    j = 1 - (i*0.01)
    k = 0.99 - (i*0.01)
    
    # Upper bound of percentile
    ub = as.numeric(unname(spatstat.geom::weighted.quantile(m.df$Input, 
                                                       m.df$Weight, 
                                                       j, 
                                                       na.rm = TRUE)))
    
    # Lower bound of percentile 
    lb = as.numeric(unname(spatstat.geom::weighted.quantile(m.df$Input, 
                                                       m.df$Weight, 
                                                       k, 
                                                       na.rm = TRUE)))
    
    # Record percentile value in percentile raster
    PercentileRaster[RasterToClassify <= ub & RasterToClassify > lb] <- j
    setTxtProgressBar(pb, i) # Update progress bar
  }
  
  PercentileRaster[is.na(PercentileRaster)] <- 0
  PercentileRaster[is.na(MaskRaster) | MaskRaster != 1] <- NA # mask classified by mask raster
  plot(PercentileRaster)
  
  return(PercentileRaster)
}