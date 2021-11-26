# Name: ht_breaks2.0 
# Description: Calculates breaks in a 1-dimensional distribution using the Head/tail breaks classification scheme 

ht_breaks <- function(x, tsh){
  # Function arguments:
  # x: Vector to classify 
  # tsh: 'Head' proportion threshold to end recursive function
  
  allheads <- c()
  
  # Breaks identification
  ht_inner <- function(x, mu){
    n <- length(x)
    mu <- c(mu, mean(x))
    
    # Clip to just 'head' values
    h <- x[x >= mean(x)]
    
    # Threshold check
    headfrac <- length(h)/n
    allheads <- c(allheads, headfrac)
    message(mean(allheads))
    
    # Recursive function 
    if(length(h) > 1 && mean(allheads) <= tsh){
      ht_inner(h, mu)
    } else mu
  }
  
  ht_inner(x, NULL)
}