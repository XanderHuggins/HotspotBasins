# RasterGridBuffer ----
# A custom buffering fuction for rasters to fill holes

RasterGridBuffer <- function(RastertoBuffer){
  # @RastertoBuffer: Rater to buffer by 1 grid cell, using logic shown below
  # Created to enlarge the JPL land mask 
  
  RastertoBuffer[is.na(RastertoBuffer[])] <- 0
  RastertoBuffer[RastertoBuffer[] != 0] <- 1
  
  Input.as.matrix <- as.matrix(RastertoBuffer)
  
  # Buffer JPL land mask by 1 grid cell to make inclusive of other data products 
  Raster.buffer <- matrix(nrow = nrow(Input.as.matrix),
                          ncol = ncol(Input.as.matrix))
  
  for (i in 2:(nrow(as.matrix(Input.as.matrix))-1)) {
    for (j in 2:(ncol(as.matrix(Input.as.matrix))-1)) {
      surrounding_sum = 0
      surrounding_sum <- sum(Input.as.matrix[i+1,j],  
                             Input.as.matrix[i+1,j+1], 
                             Input.as.matrix[i,j+1],
                             Input.as.matrix[i-1,j+1],
                             Input.as.matrix[i-1,j],   
                             Input.as.matrix[i-1,j-1],
                             Input.as.matrix[i,j-1],  
                             Input.as.matrix[i+1,j-1], na.rm = T)
      
      
      if(surrounding_sum != 0 & Input.as.matrix[i,j] == 0){
        Raster.buffer[i,j] <- 1
      } else {
        Raster.buffer[i,j] <- Input.as.matrix[i,j]
      }
    }
  }
  Raster.buffer <- raster(Raster.buffer)
  crs(Raster.buffer) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ 
                           towgs84=0,0,0"
  extent(Raster.buffer) <- extent(RastertoBuffer)
  return(Raster.buffer)
}