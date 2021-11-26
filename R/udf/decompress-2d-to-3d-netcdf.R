# Name: Decompress2dto3dNCDf 
# Description: Converts a 2-dimensional netCDF to 3-dimensional netCDF and exports as a GeoTIFF

Decompress2dto3dNCDf <- function(Raw2dNCDF, varname, var.unit, CRS.set, WriteLocation) {
  # Function arguments:
  # Raw2dNCDF: The 2-dimensional netCDF 
  # varname: Variable name within netCDF to extract
  # var.unit: Units of varname
  # CRS.set: Coordinate reference system to use
  # WriteLocation: Path to write 3-dimensional netCDF 
  
  # Import netCDF and pull attributes
  TempNCDF <- Raw2dNCDF
  lon <- ncdf4::ncvar_get(TempNCDF, "lon")
  lat <- ncdf4::ncvar_get(TempNCDF, "lat")
  vname <- varname
  message("Data attributes found")
  
  # Pull data specified from netCDF
  DataExtr <- ncvar_get(TempNCDF, vname)
  message("Data extracted")
  
  # Create empty raster stack to populate
  target_stack = stack()
  
  # Loop through each layer (e.g. month) of the dataset, converting 2d netCDF to 3d netCDF and export
  for (j in 1:dim(DataExtr)[2]) { 
    Template.matrix <- matrix(nrow = 360, ncol = 720)
    Template.matrix[] <- NA
    
    for (i in 1:dim(DataExtr)[1]) { 
      x = (lon[i] + 180.25)/0.5
      y = (lat[i] + 90.25)/0.5
      d = DataExtr[i, j]
      Template.matrix[y, x] <- d
    }
    
    Filled.ras <- raster(Template.matrix)
    Filled.ras <- flip(Filled.ras, direction = 'y')
    extent(Filled.ras) <- c(-180, 180, -90, 90)
    crs(Filled.ras) <- CRS("+init=epsg:4326")
    target_stack = stack(target_stack, Filled.ras)
    
    if (j %% 20 == 0) { print(j)} # Provides progress update
    
  }
  message("Data decompressed")
  
  # Name each layer and set CRS
  names(target_stack) <- sprintf(paste0(varname,"_%d", sep = ""), seq(1:dim(DataExtr)[2]))
  crs(target_stack) <- CRS.set
  
  # Write raster stack
  writeRaster(target_stack, WriteLocation, 
              overwrite=TRUE, 
              format="CDF", 
              varname= vname, 
              varunit= var.unit, 
              xname="lon", 
              yname="lat")
  
  message("3d netcdf stack written")
}