# Decompress2dto3dNCDf ----
# Convert 2-dimensional netCDF to 3-dimensional netCDF and export as raster

Decompress2dto3dNCDf <- function(Raw2dNCDF, varname, var.unit, CRS.set, WriteLocation) {
  # @Raw2dNCDF: The raw 2-dimensional netcdf 
  # @varname: Variable name of netcdf to extract
  # @var.unit: Variable units
  # @CRS.set: Coordinate reference system to use
  # @WriteLocation: Path to write 3-dimensional netcdf.
  
  # load netcdf and attributes
  TempNCDF <- Raw2dNCDF
  lon <- ncdf4::ncvar_get(TempNCDF, "lon")
  lat <- ncdf4::ncvar_get(TempNCDF, "lat")
  vname <- varname
  message("Data attributes found")
  
  # extract data
  DataExtr <- ncvar_get(TempNCDF, vname)
  message("Data extracted")
  
  # create empty raster stack to populate
  target_stack = stack()
  
  # loop through each month, converting 2d ncdf to 3d
  for (j in 1:dim(DataExtr)[2]) { # repeat for 480 months in dataset
    Template.matrix <- matrix(nrow = 360, ncol = 720)
    Template.matrix[] <- NA
    
    for (i in 1:dim(DataExtr)[1]) { # loop through each entry in 2d ncdf
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
    
    if (j %% 20 == 0) { print(j)}
    
  }
  message("Data decompressed")
  
  names(target_stack) <- sprintf(paste0(varname,"_%d", sep = ""), seq(1:dim(DataExtr)[2]))
  crs(target_stack) <- CRS.set
  
  writeRaster(target_stack, WriteLocation, 
              overwrite=TRUE, 
              format="CDF", 
              varname= vname, 
              varunit= var.unit, 
              xname="lon", 
              yname="lat")
  
  message("3d netcdf stack written")
}