
## Detect density peaks in allele fraction distributions

setwd("/Users/praveenp/Desktop/Chapter_2_analysis/27_10_25/Auto_allo_allele_freq/")

library(tidyverse)

# Read processed allele-fraction data
ad_long <- readRDS("allele_fraction_data.rds")

# Detect density maxima for each ploidy
dens_peaks <- ad_long %>%
  group_by(ploidy) %>%
  do({
    d <- density(
      .$AF,
      na.rm = TRUE,
      from = 0,
      to = 1
    )

    peaks_idx <- which(
      diff(sign(diff(d$y))) == -2
    ) + 1

    data.frame(
      x = d$x[peaks_idx],
      y = d$y[peaks_idx]
    )
  })

print(dens_peaks)

# Save detected peaks
write.csv(
  dens_peaks,
  "allele_fraction_density_peaks.csv",
  row.names = FALSE
)
