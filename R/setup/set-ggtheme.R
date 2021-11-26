# Set additional arguments for ggplot's theme_minimal()
theme1 = theme_minimal()+
  theme(legend.title = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        axis.line = element_line(size = 1),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())