# StressDry; UnstressDry; StressWet; UnstressWet ----
# Functions to calculate social-ecological activity statistics for various FW stress and TWS trend combinations

StressDry <- function(df.in, col.in, norm.by) {
  glob.sum <- df.in %>% pull(col.in) %>% sum(na.rm = T)
  
  tot <- df.in %>% filter(tws < -3 & fwstrs > 0.1) %>% pull(col.in) %>% 
    sum(na.rm = T)/norm.by
  
  pct <- 100 * (tot*norm.by/glob.sum)
  
  message(paste0("Total is: ", round(tot, 2), "by", norm.by, sep = ""))
  message(paste0("Pct is: ", round(pct, 1), sep = ""))
  
}

UnstressDry <- function(df.in, col.in, norm.by) {
  glob.sum <- df.in %>% pull(col.in) %>% sum(na.rm = T)
  
  tot <- df.in %>% filter(tws < -3 & fwstrs <= 0.1) %>% pull(col.in) %>% 
    sum(na.rm = T)/norm.by
  
  pct <- 100 * (tot*norm.by/glob.sum)
  
  message(paste0("Total is: ", round(tot, 2), "by", norm.by, sep = ""))
  message(paste0("Pct is: ", round(pct, 1), sep = ""))
  
}

StressWet <- function(df.in, col.in, norm.by) {
  glob.sum <- df.in %>% pull(col.in) %>% sum(na.rm = T)
  
  tot <- df.in %>% filter(tws > 3 & fwstrs > 0.1) %>% pull(col.in) %>% 
    sum(na.rm = T)/norm.by
  
  pct <- 100 * (tot*norm.by/glob.sum)
  
  message(paste0("Total is: ", round(tot, 2), "by", norm.by, sep = ""))
  message(paste0("Pct is: ", round(pct, 1), sep = ""))
  
}

UnstressWet <- function(df.in, col.in, norm.by) {
  glob.sum <- df.in %>% pull(col.in) %>% sum(na.rm = T)
  
  tot <- df.in %>% filter(tws > 3 & fwstrs <= 0.1) %>% pull(col.in) %>% 
    sum(na.rm = T)/norm.by
  
  pct <- 100 * (tot*norm.by/glob.sum)
  
  message(paste0("Total is: ", round(tot, 2), "by", norm.by, sep = ""))
  message(paste0("Pct is: ", round(pct, 1), sep = ""))
  
}