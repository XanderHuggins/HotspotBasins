# ModerateVul; HighVul; VeryHighVul ----
# Functions to calculate social-ecological activity statistics for vulnerability classes 

ModerateVul <- function(df.in, col.in, norm.by) {
  glob.sum <- df.in %>% pull(col.in) %>% sum(na.rm = T)
  
  tot <- df.in %>% filter(overallprod >= htb_o[1] & overallprod < htb_o[2]) %>% 
    pull(col.in) %>% sum(na.rm = T)/norm.by
  
  pct <- 100 * (tot*norm.by/glob.sum)
  
  message(paste0("Total is: ", round(tot, 2), "by", norm.by, sep = ""))
  message(paste0("Pct is: ", round(pct, 1), sep = ""))
  
}

HighVul <- function(df.in, col.in, norm.by) {
  glob.sum <- df.in %>% pull(col.in) %>% sum(na.rm = T)
  
  tot <- df.in %>% filter(overallprod >= htb_o[2] & overallprod < htb_o[3]) %>% 
    pull(col.in) %>% sum(na.rm = T)/norm.by
  
  pct <- 100 * (tot*norm.by/glob.sum)
  
  message(paste0("Total is: ", round(tot, 2), "by", norm.by, sep = ""))
  message(paste0("Pct is: ", round(pct, 1), sep = ""))
  
}

VeryHighVul <- function(df.in, col.in, norm.by) {
  glob.sum <- df.in %>% pull(col.in) %>% sum(na.rm = T)
  
  tot <- df.in %>% filter(overallprod >= htb_o[3]) %>% 
    pull(col.in) %>% sum(na.rm = T)/norm.by
  
  pct <- 100 * (tot*norm.by/glob.sum)
  
  message(paste0("Total is: ", round(tot, 2), "by", norm.by, sep = ""))
  message(paste0("Pct is: ", round(pct, 1), sep = ""))
  
}