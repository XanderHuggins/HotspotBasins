# tmap_clipproj ----
# clip and reproject a raster for plotting in tmap in Robinson projection

tmap_clipproj <- function(InputRaster){
  # @InputRaster: Raster to convert
  
  clip.r <- raster::crop(InputRaster, extent(-179, 180, -60, 88))
  rpj.r = projectRaster(clip.r, crs = crs("+proj=robin"), method = 'ngb')
  return(rpj.r)
}