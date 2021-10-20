# ht_breaks2.0 ----
# Head/tail breaks classification scheme 

ht_breaks2.0 <- function(x, tsh){
  # @x: vector to classify
  # @tsh: 'head' proportion threshold to end recursive function
  
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