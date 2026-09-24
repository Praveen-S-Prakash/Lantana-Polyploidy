
## Calculate allele fractions from AD values in a mixed-ploidy VCF

setwd("/Users/praveenp/Desktop/Chapter_2_analysis/27_10_25/Auto_allo_allele_freq/")

library(vcfR)
library(tidyverse)
library(data.table)
library(ggplot2)
library(scales)

# Settings
vcf_file <- "ddRAD_freebayse_mixed_ploidy_noIndel_biallelic_depth_GQ30_Q30_mac3_snp_maxM98.vcf"

minDP <- 15
maxDP <- 500

ploidy_assign <- read.csv(
  "Smple_ploidy_in_vcf_order.csv",
  stringsAsFactors = FALSE
)

# Read VCF
vcf <- read.vcfR(vcf_file, verbose = FALSE)

# Keep only biallelic SNPs
is_bial <- is.biallelic(vcf)
vcf <- vcf[is_bial, ]

# Extract allele depths
ad_mat <- extract.gt(
  vcf,
  element = "AD",
  as.numeric = FALSE
)

# Convert AD matrix to long format
variants <- rownames(ad_mat)
samples <- colnames(ad_mat)

ad_long <- data.table(
  expand.grid(
    variant = variants,
    sample = samples,
    stringsAsFactors = FALSE
  )
)

ad_values <- as.vector(ad_mat)

ad_long[, AD := ad_values]

rm(ad_values)
gc()

# Calculate allele fraction
ad_long[, AF := vapply(
  AD,
  FUN = function(x) {
    if (is.na(x) || x == "./." || x == ".") {
      return(NA_real_)
    }

    sp <- strsplit(x, ",")[[1]]

    if (length(sp) < 2) {
      return(NA_real_)
    }

    r <- as.numeric(sp[1])
    a <- as.numeric(sp[2])

    if (is.na(r) || is.na(a)) {
      return(NA_real_)
    }

    denom <- r + a

    if (denom == 0) {
      return(NA_real_)
    } else {
      return(a / denom)
    }
  },
  FUN.VALUE = 0.0
)]

# Attach ploidy information
ploidy_map <- setNames(
  ploidy_assign$ploidy,
  ploidy_assign$sample
)

ad_long[, ploidy := ploidy_map[sample]]

# Calculate depth from AD
ad_long[, DP := vapply(
  AD,
  FUN = function(x) {
    if (is.na(x) || x == "./." || x == ".") {
      return(NA_real_)
    }

    sp <- strsplit(x, ",")[[1]]

    if (length(sp) < 2) {
      return(NA_real_)
    }

    r <- as.numeric(sp[1])
    a <- as.numeric(sp[2])

    if (is.na(r) || is.na(a)) {
      return(NA_real_)
    } else {
      return(r + a)
    }
  },
  FUN.VALUE = 0.0
)]

# Filter by depth
ad_long <- ad_long[
  DP >= minDP &
    (is.na(maxDP) | DP <= maxDP)
]

# Keep heterozygous-like sites
eps <- 1e-6

ad_long <- ad_long[
  AF > eps &
    AF < (1 - eps)
]

# Save processed data
saveRDS(
  ad_long,
  file = "allele_fraction_data.rds"
)
