# BinMaker ----
# Jenks natural breaks

BinMaker <- function(VectorSet, Threshold) {
  # @VectorSet: Vector to identify natural breaks of
  # @Threshold: Goodness of fit variance threshold to stop iteration.
  
  # can't have less than 2 bins
  n_bin = 2
  
  repeat{
    CIobj = classInt::classIntervals(var = VectorSet,
                                     n = n_bin,
                                     style = "jenks")
    
    GoF = classInt::jenks.tests(CIobj)
    
    if(GoF[2] >= Threshold){
      break
    }
    
    n_bin = n_bin + 1
    
  }
  
  return(GoF)
}