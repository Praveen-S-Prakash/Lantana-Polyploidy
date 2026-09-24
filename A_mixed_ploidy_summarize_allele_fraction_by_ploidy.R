
## Summarize allele fractions by ploidy

setwd("/Users/praveenp/Desktop/Chapter_2_analysis/27_10_25/Auto_allo_allele_freq/")

library(tidyverse)
library(data.table)

# Read processed allele-fraction data
ad_long <- readRDS("allele_fraction_data.rds")

# Summaries per ploidy
summ <- ad_long %>%
  filter(!is.na(.data$ploidy)) %>%
  group_by(.data$ploidy) %>%
  summarise(
    N_sites = n(),
    mean_AF = mean(AF, na.rm = TRUE),
    median_AF = median(AF, na.rm = TRUE),
    sd_AF = sd(AF, na.rm = TRUE),
    prop_AF_approx_025 = mean(
      abs(AF - 0.25) < 0.05,
      na.rm = TRUE
    ),
    prop_AF_approx_05 = mean(
      abs(AF - 0.50) < 0.05,
      na.rm = TRUE
    ),
    prop_AF_approx_075 = mean(
      abs(AF - 0.75) < 0.05,
      na.rm = TRUE
    )
  )

print(summ)

# Save summary
write.csv(
  summ,
  "allele_fraction_summary_by_ploidy.csv",
  row.names = FALSE
)
