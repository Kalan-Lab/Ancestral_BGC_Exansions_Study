library(ggplot2)

dat <- read.table("Order_Codoff_Percentiles.txt", header=T, sep='\t')

pdf("Length_vs_Codoff.pdf", height=5, width=4)

# 2D density map
ggplot(dat, aes(x=Length, y=Codoff_Discordance_Percentile)) +
  geom_bin_2d() +
  scale_y_log10() +
  scale_x_log10() +
  labs(title = "",
       x = "BGC Length",
       y = "Codoff Discordance Percentile") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_fill_viridis_c()

dev.off()
