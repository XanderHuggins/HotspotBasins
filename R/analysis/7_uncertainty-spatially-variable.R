################################################################################
# @Manuscript - "Hotspots of social and ecological impacts from freshwater stress and storage loss" (Huggins et al.) 
# @Description - Perform spatially variable uncertainty analysis.
################################################################################

# Load 'here' library for easy path management.  
library(here)

# Import all setup and user-defined functions in R/setup and R/udfs folders
invisible(sapply(paste0(here("R/setup"), "/", list.files(here("R/setup"))), source)) 
invisible(sapply(paste0(here("R/udfs"), "/", list.files(here("R/udfs"))), source))

# Set a seed for reproducability
set.seed(47274497) #T9 of GSAS+GIWS

# Import data ----
cons <- raster(here('Data/Dimensions/TotalWithdrawls_2010.tif'))
roff <- raster(here('Data/GSCD_Qmean_0d5.tif'))
tws  <- raster(here('Data/Rodell_GRACE_TWSt_raw.tif')) * 10 #cm to mm
dis  <- raster(here('Data/HyBas4_cleaned.tif'))
efn <- raster(here("Data/Dimensions/EFN_0d5.tif"))
vsi <- raster(here("Data/Dimensions/VSI_aetw_0d5.tif"))
ac15 <- raster(here::here("Data", "AdaptiveCap2015.tif"))

# Calculate basin values for all input datasets considered in uncertainty analysis

# Establish constant extent
dis.ext <- dis
dis.ext[dis.ext >= 0] <- 1

# Generate percentile-based ecological sensitivity inputs
ef.ptl <- RasterAreaPercentiles(RasterToClassify = efn,
                                WeightRaster = WGS84_areaRaster(0.5),
                                MaskRaster = dis.ext,clipToExtent = "clip",
                                CRS.set = crs(dis.ext), ext.set = extent(dis.ext))
ef.ptl[ef.ptl > 0] <- (1 - ef.ptl[ef.ptl > 0]) + 0.01 # invert as high depths are less sensitive

vsi[is.na(vsi) | dis.ext != 1] <- NA
vsi.ptl <- RasterAreaPercentiles(RasterToClassify = vsi,
                                 WeightRaster = WGS84_areaRaster(0.5),
                                 MaskRaster = dis.ext,clipToExtent = "clip",
                                 CRS.set = crs(dis.ext), ext.set = extent(dis.ext))

# Invert adaptive capacity
ac <- 1 - ac15

# Vectorize rasters for faster computation
c_df <- raster::stack(dis, cons, roff, tws, vsi.ptl, ef.ptl, ac, WGS84_areaRaster(0.5)) %>% 
  as.data.frame() %>%  
  set_colnames(c('dis', 'cons', 'roff', 'tws', 'vsi', 'efn', 'ac', 'area'))
c_df <- c_df[complete.cases(c_df$dis),]

# Remove basins that do not have adaptive capacity data 
disclude <- c(1130, 1520, 1840, 2170, 2360, 2450, 2530, 3550, 4240,
              4370, 5320, 5330, 5370, 5430, 5520, 5530, 5540, 5740,
              6540, 6730, 7220, 7610, 7630, 7650, 7710, 7740)

for (i in 1:length(disclude)) {
  c_df <- c_df %>% filter(dis != disclude[i])
}

# Calculate average values per basin of all inputs 
c_df <- c_df %>% 
  group_by(dis) %>% 
  summarise(
    cons = weighted.mean(x = cons, w = area, na.rm = T),
    roff = weighted.mean(x = roff, w = area, na.rm = T),
    tws  = weighted.mean(x = tws, w = area, na.rm = T),
    efn = weighted.mean(x = efn, w = area, na.rm = T),
    vsi = weighted.mean(x = vsi, w = area, na.rm = T),
    ac   = weighted.mean(x = ac, w = area, na.rm = T),
  )

# Generate empty columns for vulnerability inputs
c_df$fws_ind <- rep(NA, nrow(c_df))
c_df$tws_ind <- rep(NA, nrow(c_df))
c_df$cind <- rep(NA, nrow(c_df))
c_df$overallsens   <- rep(NA, nrow(c_df))
c_df$overallprod   <- rep(NA, nrow(c_df))
c_df$hotsp  <- rep(NA, nrow(c_df))

# Calculate ecological sensitivity indicator 
c_df$ecosens <- (c_df$vsi + c_df$efn)/2
c_df$ecosens <- c_df$ecosens/max(c_df$ecosens, na.rm = T)


# Begin sensitivity analysis ----

# Establish number of realizations
runs = 1e4

# Initialize data frame to hold sensitivity results
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
bt_df <- matrix(nrow = runs, ncol = nrow(c_df)) %>% as.data.frame() %>% 
  set_colnames(c(sprintf("basin[%s]",seq(1:nrow(c_df)))))

# Initiate dataframe to track which realizations each basin is identified as a hotspot basin
bh_df <- matrix(nrow = runs, ncol = nrow(c_df)) %>% as.data.frame() %>% 
  set_colnames(c(sprintf("basin[%s]",seq(1:nrow(c_df)))))

# Set the first perturbation set to 0 for all variables to replicate the original analysis
Cptb[1] <- Qptb[1] <- Tptb[1] <- Eptb[1] <- Vptb[1] <- Aptb[1] <- 0

# Loop through all perturbation sets

for (i in 1:nrow(s_df)) {

  t_df <- c_df # Set temporary calculation dataframe
  
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
    t_df$efn[k]  <- max(min(t_df$efn[k] + Epsd[k], 1), 0) # Ensure staying within 0-1 range
    t_df$vsi[k]  <- max(min(t_df$vsi[k] + Vpsd[k], 1), 0) # Ensure staying within 0-1 range
    t_df$ac[k]   <- max(min(t_df$ac[k] + Apsd[k], 1), 0)  # Ensure staying within 0-1 range
  }
  
  # Calculate ecological sensitivity
  t_df$ecosens <- (t_df$vsi + t_df$efn)/2
  t_df$ecosens <- t_df$ecosens/max(t_df$ecosens, na.rm = T)
  
  # Calculate basin freshwater status
  for (j in 1:nrow(t_df)) {
    t_df$fws_ind[j] <- min((t_df$cons[j]/(0.4*t_df$roff[j])), 1)
    t_df$tws_ind[j] <- max(min((t_df$tws[j]/(0.4*t_df$roff[j])), 1), -1)*-1
    t_df$cind[j] <- max(min(( (t_df$fws_ind[j] + t_df$tws_ind[j])/2 ), 1), 0)
  }
  
  # Handle regions with earthquake interference of TWS trends
  eqbas <- c(4130, 4440, # Japan
             4451, 4430, 4454, 4457, 4457, 4459, 4460,
             5111, 5112, 5113, 5114, 5115, 5116, 5117, 5118, 5120, 5130) # Malay Pen.
  
  for (t in eqbas) {
    for (m in 1:nrow(t_df)) {
      
      if (t_df$dis[m] == t) {
        t_df$cind[m] <- t_df$fws_ind[m]  
      }
    }
  }
  
  # Calculate social-ecological vulnerability
  t_df$overallsens <- 1 - ((1-t_df$ecosens)*(1-t_df$ac)) # this is fuzzy sum
  t_df$overallprod <- t_df$cind * t_df$overallsens
  
  # Store perturbations applied as standard deviation of set
  s_df$Cptb[i] <- Cptb[i]
  s_df$Qptb[i] <- Qptb[i]
  s_df$Tptb[i] <- Tptb[i]
  s_df$Aptb[i] <- Aptb[i]
  s_df$Eptb[i] <- Eptb[i]
  s_df$Vptb[i] <- Vptb[i]
  
  # Identify vulnerability class breaks using Head/Tail breaks method
  htb <- ht_breaks2.0(t_df$overallprod, 0.8)
  
  # Document which basins are transitional and hotspot basins, per perturbation set
  bt_df[i,] <- t_df$overallprod > htb[1]
  bh_df[i,] <- t_df$overallprod > htb[2]
  
  # Track the number of transitional, and hotspot (high and very high vul.) basins per perturbation set
  s_df$NhotM[i] <- t_df %>% filter(overallprod >= htb[1] & overallprod < htb[2]) %>% nrow()
  s_df$NhotH[i] <- t_df %>% filter(overallprod >= htb[2] & overallprod < htb[3]) %>% nrow()
  s_df$NhotVH[i] <- t_df %>% filter(overallprod >= htb[3]) %>% nrow()
  
  print(i)
  
}

# Save results
write.csv(s_df, here('Data/montecaro_scalevariance-update.csv'))
write.csv(bt_df, here('Data/thresholdbasins_scalevariance-update.csv'))
write.csv(bh_df, here('Data/hotspotbasins_scalevariance-update.csv'))

# Plot and analyze sensitivity results ----

# Load sensitivity results (so to save from having to rerun sensitivity loop)
s_df <- read.csv(here('Data/montecaro_scalevariance-update.csv'))
bt_df <- read.csv(here('Data/thresholdbasins_scalevariance-update.csv'))
bh_df <- read.csv(here('Data/hotspotbasins_scalevariance-update.csv'))

# Select columns of interest to plot
p_df <- s_df %>% dplyr::select(NhotH, NhotVH, Cptb, Qptb, Tptb, Aptb, Eptb, Vptb)
names(p_df) <- c('NhotH', 'NhotVH', 'a) Consumption', 'b) Runoff', 'c) TWS', 'd) Adaptability', 'e) EFN', 'g) VSI')
p_df$hotspots <- p_df$NhotH + p_df$NhotVH # Calculate hotspot total
p_df <- p_df[,-c(1,2)] # Drop un-needed columns

# Generate facet plot ----
sensplot <- 
  p_df %>% 
  gather(-hotspots, key = "var", value = "value") %>%
  ggplot(aes(x = value, y = hotspots)) +
  geom_point(alpha = 0.05, size = 1) + 
  geom_point(data = NULL, aes(x = 0, y = 172), color = "red", size = 1.5) +
  # geom_smooth(alpha = 0.5) +
  facet_wrap(~ var, nrow = 1) +
  scale_y_continuous(limits = c(0, 300)) +
  scale_x_continuous(limits = c(0, 0.2)) +
  theme_bw() +
  theme(panel.spacing.x = unit(1, "lines"),
        axis.text.x = element_text(size = 6)) +
  coord_cartesian(expand = c(0,0), clip = "off")
sensplot

ggsave(filename = "C:/Users/xande/Desktop/sens_scatter_varscale.pdf", plot = sensplot,
       device = 'pdf', width = 190, height = 60, units = 'mm', dpi = 400)


# Generate plots showing ranked frequency of basins identified as transitional or hotspot basins ----

# Initiate dataframe to store frequency results
summary_df <- matrix(nrow = ncol(bt_df), ncol = 3) %>% 
  as.data.frame() %>% set_colnames(c('BasinID', 'Freq', 'Actual'))

# Count number of perturbation realizations each basin is identified as a transitional or hotspot basin
for (i in 1:nrow(summary_df)-1) {
  summary_df$BasinID[i] <- colnames(bt_df)[i+1]
  summary_df$Freq[i] <- sum(bt_df[,i+1], na.rm = T)/nrow(bt_df)
  summary_df$Actual[i] <- bt_df[1, i+1]
}

# Order data frame based on frequency and provide rank no.
summary_df <- summary_df[order(-summary_df$Freq),]
summary_df$rank <- seq(1, nrow(summary_df), 1)

# Generate plot
ggplot(summary_df, aes(x=rank, y=Freq, fill=as.factor(Actual))) + 
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = c('grey', '#B99232')) +
  coord_cartesian(xlim = c(0, 1204), ylim = c(0, 1), expand = c(0, 0), clip = "off") +
  scale_y_continuous(breaks = seq(0, 1, 0.1), 
                     limits = c(0, 1)) + 
  theme1 + theme(axis.ticks.x = element_line(size = 1)) +
  geom_hline(yintercept = 0.5, col = 'black', lwd = 0.5, linetype = "dashed") +
  geom_hline(yintercept = 0.95, col = 'black', lwd = 0.5, linetype = "dashed")

ggsave(filename = "C:/Users/xande/Desktop/threshold_frequency.pdf", 
       plot = last_plot(), device = 'pdf', width = 190/1.5, height = 100/1.5, 
       units = 'mm', dpi = 400)

# Number of transitional+hotspots identified as so in > 95% of perturbations
summary_df %>% filter(Actual == 1) %>% filter(Freq >= 0.95) %>% nrow()

# Number of transitional+hotspots identified as so in > 50% of perturbations
summary_df %>% filter(Actual == 1) %>% filter(Freq >= 0.50) %>% nrow()

# Number of transitional+hotspots identified as so in < 50% of perturbations
summary_df %>% filter(Actual == 1) %>% filter(Freq < 0.50) %>% nrow()

# Number of *non* transitional+hotspots identified as so in > 50% of perturbations
summary_df %>% filter(Actual == 0) %>% filter(Freq >= 0.50) %>% nrow()
  

# Generate plots showing ranked frequency of basins identified as hotspot basins ----

# Initiate dataframe to store frequency results
summary_df <- matrix(nrow = ncol(bh_df), ncol = 3) %>% 
  as.data.frame() %>% set_colnames(c('BasinID', 'Freq', 'Actual'))

# Count number of perturbation realizations each basin is identified as a hotspot basin
for (i in 1:nrow(summary_df)-1) {
  summary_df$BasinID[i] <- colnames(bh_df)[i+1]
  summary_df$Freq[i] <- sum(bh_df[,i+1], na.rm = T)/nrow(bh_df)
  summary_df$Actual[i] <- bh_df[1, i+1]
}

# Order data frame based on frequency and provide rank no.
summary_df <- summary_df[order(-summary_df$Freq),]
summary_df$rank <- seq(1, nrow(summary_df), 1)

# Generate plot
ggplot(summary_df, aes(x=rank, y=Freq, fill=as.factor(Actual))) + 
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = c('grey', '#7E1900')) +
  coord_cartesian(xlim = c(0, 1204), ylim = c(0, 1), expand = c(0, 0), clip = "off") +
  scale_y_continuous(breaks = seq(0, 1, 0.1), 
                     limits = c(0, 1)) + 
  theme1 + theme(axis.ticks.x = element_line(size = 1)) +
  geom_hline(yintercept = 0.5, col = 'black', lwd = 0.5, linetype = "dashed")+
  geom_hline(yintercept = 0.95, col = 'black', lwd = 0.5, linetype = "dashed")

ggsave(filename = "C:/Users/xande/Desktop/hotspot_frequency.pdf", 
       plot = last_plot(), device = 'pdf', width = 190/1.5, height = 100/1.5, 
       units = 'mm', dpi = 400)

# Number of hotspots identified as so in > 95% of perturbations
summary_df %>% filter(Actual == 1) %>% filter(Freq >= 0.95) %>% nrow()

# Number of hotspots identified as so in > 50% of perturbations
summary_df %>% filter(Actual == 1) %>% filter(Freq >= 0.50) %>% nrow()

# Number of hotspots identified as so in < 50% of perturbations
summary_df %>% filter(Actual == 1) %>% filter(Freq < 0.50) %>% nrow()

# Number of *non* hotspots identified as so in > 50% of perturbations
summary_df %>% filter(Actual == 0) %>% filter(Freq >= 0.50) %>% nrow()