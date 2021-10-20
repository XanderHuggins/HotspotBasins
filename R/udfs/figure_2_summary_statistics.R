# Modstress_badAC; Modstress_goodAC; Highstress_badAC; Highstress_goodAC ----
# Functions to calculate social-ecological activity statistics for various basin freshwater status and  adaptive capacity combinations

Modstress_badAC <- function(df.in, col.in, norm.by) {
  glob.sum <- df.in %>% pull(col.in) %>% sum(na.rm = T)
  
  tot <- df.in %>% filter(cind > ptls$ind[67] & cind <= ptls$ind[80] & ac > ptls$ac[67]) %>% 
    pull(col.in) %>% sum(na.rm = T)/norm.by
  
  pct <- 100 * (tot*norm.by/glob.sum)
  
  message(paste0("Total is: ", round(tot, 2), "by", norm.by, sep = ""))
  message(paste0("Pct is: ", round(pct, 1), sep = ""))
  
}

Modstress_goodAC <- function(df.in, col.in, norm.by) {
  glob.sum <- df.in %>% pull(col.in) %>% sum(na.rm = T)
  
  tot <- df.in %>% filter(cind > ptls$ind[67] & cind <= ptls$ind[80] & ac <= ptls$ac[67]) %>% 
    pull(col.in) %>% sum(na.rm = T)/norm.by
  
  pct <- 100 * (tot*norm.by/glob.sum)
  
  message(paste0("Total is: ", round(tot, 2), "by", norm.by, sep = ""))
  message(paste0("Pct is: ", round(pct, 1), sep = ""))
  
}

Highstress_badAC <- function(df.in, col.in, norm.by) {
  glob.sum <- df.in %>% pull(col.in) %>% sum(na.rm = T)
  
  tot <- df.in %>% filter(cind > ptls$ind[80] & ac > ptls$ac[67]) %>% 
    pull(col.in) %>% sum(na.rm = T)/norm.by
  
  pct <- 100 * (tot*norm.by/glob.sum)
  
  message(paste0("Total is: ", round(tot, 2), "by", norm.by, sep = ""))
  message(paste0("Pct is: ", round(pct, 1), sep = ""))
  
}

Highstress_goodAC <- function(df.in, col.in, norm.by) {
  glob.sum <- df.in %>% pull(col.in) %>% sum(na.rm = T)
  
  tot <- df.in %>% filter(cind > ptls$ind[80] & ac <= ptls$ac[67]) %>% 
    pull(col.in) %>% sum(na.rm = T)/norm.by
  
  pct <- 100 * (tot*norm.by/glob.sum)
  
  message(paste0("Total is: ", round(tot, 2), "by", norm.by, sep = ""))
  message(paste0("Pct is: ", round(pct, 1), sep = ""))
  
}
