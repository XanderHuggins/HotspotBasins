# Name: figure-S8-plotting+stats.R
# Description: Plot scatter plot of uncertainty results, ordered bar graphs of hotspot identification frequency of individual basins, and calculate associated summary statistics

# Use 'here' package for easy path management
library(here)

# Import all setup and user-defined functions in R/setup and R/udfs folders
invisible(sapply(paste0(here("R/setup"), "/", list.files(here("R/setup"))), source)) 
invisible(sapply(paste0(here("R/udf"), "/", list.files(here("R/udf"))), source))

# Import coastlines and basin shapefile with plotting data
s_df <- read.csv(here('Data/mc_basin_all_uniform.csv'))
bt_df <- read.csv(here('Data/mc_basin_threshold_uniform.csv'))
bh_df <- read.csv(here('Data/mc_basin_hotspot_uniform.csv'))

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
  geom_point(data = NULL, aes(x = 0, y = 168), color = "red", size = 1.5) +
  # geom_smooth(alpha = 0.5) +
  facet_wrap(~ var, nrow = 1) +
  scale_y_continuous(limits = c(0, 300)) +
  scale_x_continuous(limits = c(-0.2, 0.2)) +
  theme_bw() +
  theme(panel.spacing.x = unit(1, "lines"),
        axis.text.x = element_text(size = 6)) +
  coord_cartesian(expand = c(0,0))
sensplot

ggsave(filename = here::here("plot_outputs/FigS8_sens_uniform_scatter.pdf"), plot = sensplot,
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

ggsave(filename = here::here("plot_outputs/FigS8_frequency_uniform_perturbations.pdf"), 
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

ggsave(filename = here::here("plot_outputs/FigS8_frequency_uniform_perturbation.pdf"), 
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
