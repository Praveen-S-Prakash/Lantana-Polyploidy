
## Plot allele fraction distributions by ploidy

setwd("/Users/praveenp/Desktop/Chapter_2_analysis/27_10_25/Auto_allo_allele_freq/")

library(ggplot2)

# Read processed allele-fraction data
ad_long <- readRDS("allele_fraction_data.rds")

ploidy_colors <- c(
  "2" = "#8863ce",
  "3" = "#edf41b",
  "4" = "#33d5db",
  "6" = "#ec501b"
)

ad_long$ploidy <- factor(
  ad_long$ploidy,
  levels = sort(unique(ad_long$ploidy))
)

# Faceted histogram and density
p <- ggplot(
  ad_long,
  aes(x = AF, color = ploidy, fill = ploidy)
) +
  geom_histogram(
    position = "identity",
    alpha = 0.35,
    bins = 80
  ) +
  geom_density(
    alpha = 0.0,
    adjust = 1.0,
    size = 0.8
  ) +
  facet_wrap(
    ~ploidy,
    ncol = 1,
    scales = "free_y"
  ) +
  scale_color_manual(values = ploidy_colors) +
  scale_fill_manual(values = ploidy_colors) +
  labs(
    x = "ALT allele fraction (alt/(ref+alt))",
    y = "Count of genotype calls",
    title = "Allele balance by ploidy"
  ) +
  theme_minimal(base_size = 14)

print(p)

ggsave(
  "allele_fraction_by_ploidy_facet_coloured.pdf",
  plot = p,
  width = 6,
  height = 10
)

# Combined density plot
p2 <- ggplot(
  ad_long,
  aes(x = AF, color = ploidy)
) +
  geom_density(size = 1) +
  scale_color_manual(values = ploidy_colors) +
  labs(
    x = "ALT allele fraction",
    y = "Density",
    title = "Overlay allele-frequency density by ploidy"
  ) +
  theme_minimal(base_size = 14)

print(p2)

ggsave(
  "allele_fraction_by_ploidy_overlay_coloured.pdf",
  plot = p2,
  width = 7,
  height = 5
)
