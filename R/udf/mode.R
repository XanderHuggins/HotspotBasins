# Name: getmode 
# Description: Calculates mode of provided vector

getmode <- function(v) {
  # Function arguments:
  # v: Vector to calculate mode of
  
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}