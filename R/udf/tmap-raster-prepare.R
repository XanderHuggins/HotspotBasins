# Name: tmap_clipproj 
# Description: Prepares a global raster in WGS84 for plotting using tmap by cropping and reprojecting to the Robinson projection   

tmap_clipproj <- function(InputRaster){
  # Function arguments:
  # InputRaster: Raster to convert
  
  clip.r <- raster::crop(InputRaster, extent(-179, 180, -60, 88))
  rpj.r = projectRaster(clip.r, crs = crs("+proj=robin"), method = 'ngb')
  return(rpj.r)
}