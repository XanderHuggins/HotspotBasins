# Name: 4-uncertainty-spatially-variable.R
# Description: Perform uncertainty analysis to consider the impacts of potential spatially variable uncertainty in all input datasets   

# Use 'here' package for easy path management
library(here)

# Import all setup and user-defined functions in R/setup and R/udfs folders
invisible(sapply(paste0(here("R/setup"), "/", list.files(here("R/setup"))), source)) 
invisible(sapply(paste0(here("R/udf"), "/", list.files(here("R/udf"))), source))

# Import original data
basins_data <- sf::read_sf(here('R/plotting+stats/Basin_data.shp'))
basins_data <- basins_data[complete.cases(basins_data$ss_vcls), ] 
basins_data <- basins_data %>% 
  as.data.frame() %>%  
  dplyr::select(pfaf_id, tws, wuse, roff, efn, vsi, ac) 
colnames(basins_data) <- c('pfaf_id', 'tws', 'cons', 'roff', 'efn', 'vsi', 'ac')

# Set seed for reproducibility
set.seed(47274497) #T9 of GSAS+GIWS

# Begin sensitivity analysis ----

# Establish number of realizations
runs = 1e4

# Initialize data frame to hold sensitivity results, and create perturbation sets from 
s_df <- data.frame(id = seq(1, runs, by = 1),
                   Cpsd = rep(NA, runs), 
                   Qpsd = rep(NA, runs),
                   Tpsd = rep(NA, runs),
                   Apsd = rep(NA, runs),
                   Epsd = rep(NA, runs),
                   Vpsd = rep(NA, runs),
                   NhotM = rep(NA, runs), 
                   NhotH = rep(NA, runs), 
                   NhotVH = rep(NA, runs))

# Randomly sample from uniform distribution to change variance of spatially variable perturbations to apply to each variable  
Cptb <- runif(runs, min = 0, max = 0.2)
Qptb <- runif(runs, min = 0, max = 0.2)
Tptb <- runif(runs, min = 0, max = 0.2)
Eptb <- runif(runs, min = 0, max = 0.2)
Vptb <- runif(runs, min = 0, max = 0.2)
Aptb <- runif(runs, min = 0, max = 0.2)

# Initiate dataframe to track which realizations each basin is identified as a transitional or hotspot basin
bt_df <- matrix(nrow = runs, ncol = nrow(basins_data)) %>% as.data.frame() %>% 
  set_colnames(c(sprintf("basin[%s]",seq(1:nrow(basins_data)))))

# Initiate dataframe to track which realizations each basin is identified as a hotspot basin
bh_df <- matrix(nrow = runs, ncol = nrow(basins_data)) %>% as.data.frame() %>% 
  set_colnames(c(sprintf("basin[%s]",seq(1:nrow(basins_data)))))

# Set the first perturbation set to 0 for all variables to replicate the original analysis
Cptb[1] <- Qptb[1] <- Tptb[1] <- Eptb[1] <- Vptb[1] <- Aptb[1] <- 0

# Loop through all perturbation sets

for (i in 1:nrow(s_df)) {

  t_df <- basins_data # Set temporary calculation dataframe
  
  # Generate spatially variable perturbations to apply to each basin per run, based on random sampling from normal distribution with mean = 0 and standard deviation based on multiplier randomly sampled above
  Cpsd <- rnorm(nrow(t_df), mean = 0, sd = Cptb[i]*sqrt(var(t_df$cons)))
  Qpsd <- rnorm(nrow(t_df), mean = 0, sd = Qptb[i]*sqrt(var(t_df$roff)))
  Tpsd <- rnorm(nrow(t_df), mean = 0, sd = Tptb[i]*sqrt(var(t_df$tws)))
  Epsd <- rnorm(nrow(t_df), mean = 0, sd = Eptb[i]*sqrt(var(t_df$efn)))
  Vpsd <- rnorm(nrow(t_df), mean = 0, sd = Vptb[i]*sqrt(var(t_df$vsi)))
  Apsd <- rnorm(nrow(t_df), mean = 0, sd = Aptb[i]*sqrt(var(t_df$ac)))
  
  
  # Apply perturbations to basins 
  for (k in 1:nrow(t_df)) {
    t_df$cons[k] <- max(t_df$cons[k] + Cpsd[k], 0) # Ensure no negative withdrawal rates
    t_df$roff[k] <- max(t_df$roff[k] + Qpsd[k], 0.1) # Ensure there is are no 0 streamflow values
    t_df$tws[k]  <- t_df$tws[k]  + Tpsd[k]
    t_df$efn[k]  <- max(min(t_df$efn[k] + Epsd[k], 1), 0) # Ensure within 0-1 range
    t_df$vsi[k]  <- max(min(t_df$vsi[k] + Vpsd[k], 1), 0) # Ensure within 0-1 range
    t_df$ac[k]   <- max(min(t_df$ac[k] + Apsd[k], 1), 0)  # Ensure within 0-1 range
  }
  
  # Calculate ecological sensitivity
  t_df$ecol_sens <- (t_df$vsi + t_df$efn)/2
  t_df$ecol_sens <- t_df$ecol_sens/max(t_df$ecol_sens, na.rm = T)
  
  # Calculate social sensitivity
  t_df$socl_sens <- 1 - t_df$ac
  
  # Calculate basin freshwater status
  for (j in 1:nrow(t_df)) {
    t_df$fws_ind[j] <- min((t_df$cons[j]/(0.4*t_df$roff[j])), 1)
    t_df$tws_ind[j] <- max(min((t_df$tws[j]/(0.4*t_df$roff[j])), 1), -1) * -1
    t_df$fw_status[j] <- max(min(( (t_df$fws_ind[j] + t_df$tws_ind[j])/2 ), 1), 0)
  }
  
  # Handle regions with earthquake interference of TWS trends
  eq_intfr <- readr::read_csv(file = here::here("Data/basin_lists_1.csv"),
                              col_select = 'eq_interference',
                              col_names = T, 
                              show_col_types = FALSE)
  
  for (t in 1:nrow(t_df)) { 
    if (t_df$pfaf_id[t] %in% eq_intfr$eq_interference) {
      t_df$fw_status[t] <- t_df$fws_ind[t]
    }
  }
  
  # Calculate social-ecological vulnerability
  t_df$ses_sens <- 1 - ((1-t_df$ecol_sens)*(1-t_df$socl_sens))
  t_df$ses_vuln <- t_df$fw_status * t_df$ses_sens
  
  # Store perturbations applied as standard deviation of set
  s_df$Cpsd[i] <- Cptb[i]
  s_df$Qpsd[i] <- Qptb[i]
  s_df$Tpsd[i] <- Tptb[i]
  s_df$Apsd[i] <- Aptb[i]
  s_df$Epsd[i] <- Eptb[i]
  s_df$Vpsd[i] <- Vptb[i]
  
  # Identify vulnerability class breaks using Head/Tail breaks method
  htb <- ht_breaks(t_df$ses_vuln, 0.8)
  
  # Document which basins are transitional and hotspot basins, per perturbation set
  bt_df[i,] <- t_df$ses_vuln > htb[1]
  bh_df[i,] <- t_df$ses_vuln > htb[2]
  
  # Track the number of transitional, and hotspot (high and very high vul.) basins per perturbation set
  s_df$NhotM[i] <- t_df %>% dplyr::filter(ses_vuln >= htb[1] & ses_vuln < htb[2]) %>% nrow()
  s_df$NhotH[i] <- t_df %>% filter(ses_vuln >= htb[2] & ses_vuln < htb[3]) %>% nrow()
  s_df$NhotVH[i] <- t_df %>% filter(ses_vuln >= htb[3]) %>% nrow()
  
  print(i/runs)
  
}

# Save results
write.csv(s_df, here('Data/mc_basin_all_variable.csv'))
write.csv(bt_df, here('Data/mc_basin_threshold_variable.csv'))
write.csv(bh_df, here('Data/mc_basin_hotspot_variable.csv'))