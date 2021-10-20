# getmode ----
# pulls most frequent value from vector

getmode <- function(v) {
  # @v: Vector to calculate mode of
  
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


# getmode2 ----
# pulls second most frequent value from vector 

getmode2 <- function(v) {
  uniqv <- unique(v)
  modev <- uniqv[which.max(tabulate(match(v, uniqv)))]
  
  dropmode <- v[v != modev]
  
  if (length(dropmode) > 0) {
    uniqv2 <- unique(dropmode)
    uniqv2[which.max(tabulate(match(dropmode, uniqv2)))]
  } else {
    return(NA)
  }
  
}