
## Convert mixed-ploidy VCF to diploid genotypes for STRUCTURE/ADMIXTURE

library(vcfR)
library(tidyverse)

setwd("/Users/praveenp/Desktop/Chapter_2_analysis/27_10_25/Structure_analysis/")

vcf <- read.vcfR(
  "ddRAD_freebayse_mixed_ploidy_noIndel_biallelic_depth_GQ30_Q30_mac3_snp_maxM98_Structure.vcf"
)

gt <- extract.gt(vcf, element = "GT", as.numeric = FALSE)

dip_gt <- gt

sample_two <- function(gt_string) {
  if (is.na(gt_string) || gt_string == "." || gt_string == "./.") {
    return("./.")
  }

  alleles <- unlist(strsplit(gsub("[|/]", " ", gt_string), " "))
  alleles <- alleles[alleles != ""]

  if (length(alleles) < 2) {
    return("./.")
  }

  # Sample two alleles at random
  picked <- sample(alleles, size = 2, replace = TRUE)

  # Sort so that 0/1 is not written as 1/0
  picked <- sort(picked)

  paste0(picked[1], "/", picked[2])
}

for (i in 1:ncol(gt)) {
  for (j in 1:nrow(gt)) {
    dip_gt[j, i] <- sample_two(gt[j, i])
  }
}

# Write diploidized genotypes back into the VCF
vcf@gt[, -1] <- dip_gt

write.vcf(
  vcf,
  file = "diploidized_for_admixture.vcf"
)
