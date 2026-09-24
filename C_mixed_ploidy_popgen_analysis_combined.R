## Mixed ploidy analysis Paper 2 ##

library(scales)
library(adegenet)
library(adegraphics)
library(vcfR)
library(pegas)
library(StAMPP)
library(gplots)
library(ape)

setwd("/Users/praveenp/Desktop/Chapter_2_analysis/27_10_25/")
source("./adgenet_functions_mixed.R")
##############################
### IMPORT SNP data from VCF 
##############################

vcf <- read.vcfR("ddRAD_freebayse_mixed_ploidy_noIndel_biallelic_depth_GQ30_Q30_mac3_snp_maxM98.vcf")      #,nrows=10000 # SNPs, max 10% of missing data per site (AN > 850): 1,350,328 4-fold degenerated SNPs
head(vcf)               #check
vcf@fix[1:10,1:5]      #check

#import cvs with sample and ploidy in the order
ploidy_table <- read.csv("Sample_ploidy_in_vcf_order.csv", header = TRUE)

### convert to genlight
lc.genlight <- vcfR2genlight.mixed(vcf, ploidy_table = ploidy_table)
locNames(lc.genlight) <- paste(vcf@fix[,1],vcf@fix[,2],sep="_")   # add real SNP.names
pop(lc.genlight)<-substr(indNames(lc.genlight),1,5)               # add pop names: here pop names are first 5 chars of ind name

#check the object
lc.genlight
indNames(lc.genlight)
ploidy(lc.genlight)

##############################
#  data checks & statistics
##############################

## 1) plot the entire matrix
png("glplot.png",width = 1000, height = 500, units = "px")
glPlot (lc.genlight)  # takes some time
dev.off()

### 2) N missing SNPs per sample
x <- summary(t(as.matrix(lc.genlight)))
write.table(substr(x[7,],9,100), file = "missing.persample.txt", sep = "\t", quote=F, col.names=F)  # NAs, if present, are in seventh row of summary

### 3) check persample read depth and plot it
dp <- extract.gt(vcf, element='DP', as.numeric=TRUE)
mean(colMeans(dp,na.rm=T)) #average coverage

png(file="DP.png",width=1000, height=500)
par(mar=c(8,4,1,1))
boxplot(dp, las=3, col=c("#C0C0C0", "#808080"), ylab="Read Depth (DP)", las=2, cex=0.7)
dev.off()

png(file="DP_detail.png",width=1000, height=500)
par(mar=c(8,4,1,1))
boxplot(dp, las=3, col=c("#C0C0C0", "#808080"), ylab="Read Depth (DP)", las=2, cex=0.7, ylim=c(0,500))
abline(h=4, col="red")
dev.off()

##################################
#             SFS
##################################

### 4) plot total allele frequency spectrum (AFS) of entire dataset
mySum <- glSum(lc.genlight, alleleAsUnit = TRUE) 
png("total_AFS.png",width = 1000, height = 500, units = "px")
barplot(table(mySum), col="blue", space=0, xlab="Allele counts",  main="Distribution of ALT allele counts in total dataset")
dev.off()

### 5 seperate sfs per ploidy ###########################################

#########
#Split your genlight object by ploidy
lc.di   <- lc.genlight[ploidy(lc.genlight) == 2, ]
lc.tri  <- lc.genlight[ploidy(lc.genlight) == 3, ]
lc.tetra<- lc.genlight[ploidy(lc.genlight) == 4, ]
lc.hexa <- lc.genlight[ploidy(lc.genlight) == 6, ]

ploidy_colors <- c(
  "Diploid"    = "#8863ce",  # violet
  "Triploid"   = "#edf41b",  # yellow
  "Tetraploid" = "#33d5db",  # cyan
  "Hexaploid"  = "#ec501b"   # orange-red
)

#Define a small helper function for plotting SFS
plot_SFS <- function(gl_obj, ploidy_label) {
  mySum <- glSum(gl_obj, alleleAsUnit = TRUE)
  png(paste0("AFS_", ploidy_label, ".png"), width = 1000, height = 500, units = "px")
  barplot(
    table(mySum),
    col = ploidy_colors[ploidy_label],
    border = NA,
    space = 0,
    xlab = "Allele counts",
    ylab = "Number of SNPs",
    main = paste("ALT allele frequency spectrum -", ploidy_label)
  )
  dev.off()
}

#Generate SFS for each ploidy
plot_SFS(lc.di,   "Diploid")
plot_SFS(lc.tri,  "Triploid")
plot_SFS(lc.tetra,"Tetraploid")
plot_SFS(lc.hexa, "Hexaploid")

###### 6. combined sfs of all ploidy.  ################

###############

# Here alleleic frequency is used. I dont know how we can use this

#(1) calculate the allele frequencies per SNP for each ploidy group, 
#(2) bin those frequencies, and 
#(3) plot them together with color-coded bars or density lines.

#Function to compute allele frequency (not count)
get_AF <- function(gl_obj) {
  mySum <- glSum(gl_obj, alleleAsUnit = TRUE)
  # divide by total possible alleles (ploidy × number of individuals)
  total_alleles <- sum(ploidy(gl_obj))
  af <- mySum / total_alleles
  return(af)
}

#Compute allele frequencies per ploidy group
af_di    <- get_AF(lc.di)
af_tri   <- get_AF(lc.tri)
af_tetra <- get_AF(lc.tetra)
af_hexa  <- get_AF(lc.hexa)

#Combine all into a tidy dataframe for plotting
library(ggplot2)

sfs_data <- data.frame(
  freq = c(af_di, af_tri, af_tetra, af_hexa),
  ploidy = c(
    rep("Diploid", length(af_di)),
    rep("Triploid", length(af_tri)),
    rep("Tetraploid", length(af_tetra)),
    rep("Hexaploid", length(af_hexa))
  )
)
# Define your custom colors for ploidy groups
ploidy_colors <- c(
  "Diploid"   = "#8863ce",  # purple
  "Triploid"  = "#edf41b",  # yellow
  "Tetraploid"= "#33d5db",  # teal
  "Hexaploid" = "#ec501b"   # orange
)

# Plot all SFS curves together
p_sfs <- ggplot(sfs_data, aes(x = freq, fill = ploidy, color = ploidy)) +
  geom_histogram(
    position = "identity",
    alpha = 0.5,
    bins = 50
  ) +
  scale_x_continuous(name = "Allele frequency (ALT allele)") +
  scale_y_continuous(name = "Number of SNPs") +
  scale_fill_manual(values = ploidy_colors) +
  scale_color_manual(values = ploidy_colors) +
  theme_minimal(base_size = 14) +
  ggtitle("Site Frequency Spectrum (SFS) by Ploidy Group") +
  theme(
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5)
  )

# Show and save
print(p_sfs)
ggsave("SFS_by_ploidy_coloured.pdf", plot = p_sfs, width = 8, height = 6)

##############
# 7. This is to compute some if the summary stat form the sfs data

# Combine all AFs into one data frame
sfs_data <- data.frame(
  freq = c(af_di, af_tri, af_tetra, af_hexa),
  ploidy = c(
    rep("Diploid", length(af_di)),
    rep("Triploid", length(af_tri)),
    rep("Tetraploid", length(af_tetra)),
    rep("Hexaploid", length(af_hexa))
  )
)

# Define bins for minor allele frequencies (MAF)
sfs_data$MAF <- ifelse(sfs_data$freq > 0.5, 1 - sfs_data$freq, sfs_data$freq)

# Compute summary stats per ploidy
library(dplyr)

sfs_summary <- sfs_data %>%
  group_by(ploidy) %>%
  summarise(
    N_SNPs = n(),
    Mean_AF = mean(freq, na.rm = TRUE),
    Median_AF = median(freq, na.rm = TRUE),
    SD_AF = sd(freq, na.rm = TRUE),
    Prop_Rare = mean(MAF < 0.05, na.rm = TRUE),   # % of rare variants (MAF < 0.05)
    Prop_Common = mean(MAF > 0.45, na.rm = TRUE), # % of nearly fixed/common variants
    Mean_Het = mean(2 * freq * (1 - freq), na.rm = TRUE) # expected heterozygosity
  )

print(sfs_summary)



##################################
# distance matrix for SplitsTree
##################################

# 8

# calculate matrix Nei's 1972 distance between indivs and pops
# export matrix - open it then in SplitsTree

lc.D.ind <- stamppNeisD(lc.genlight, pop = FALSE)                 # Nei's 1972 distance between indivs
stamppPhylip(lc.D.ind, file="indiv_Neis_distance_mix_ploidy.phy.dst")    

lc.D.pop <- stamppNeisD(lc.genlight, pop = TRUE)                  # Nei's 1972 distance between pops
stamppPhylip(lc.D.pop, file="pops_Neis_distance_mix_ploidy.phy.dst") 

###### 9 distance based on ploidy

# Make a copy of your genlight
lc.genlight_copy <- lc.genlight

# Set ploidy as population only in the copy
pop(lc.genlight_copy) <- ploidy(lc.genlight_copy)

# Compute Nei’s distances between ploidy groups
lc.D.ploidy <- stamppNeisD(lc.genlight_copy, pop = TRUE)
stamppPhylip(lc.D.ploidy, file="ploidy_Neis_distance.phy.dst")


##################################
#### NJ tree  ######
##################################
# 10

# NJ tree
nj.tree <- nj(lc.D.ind)   # nj() accepts 'dist' object

# simple tree
plot(nj.tree, main = "NJ Tree of Individuals", cex = 0.7)

# tip labels colored by ploidy
ploidy_factor <- as.factor(ploidy(lc.genlight))  # 2x,3x,4x,6x
tip.cols <- c("blue", "red", "green", "purple")   # assign colors to ploidy levels
tip.color <- tip.cols[as.numeric(ploidy_factor)]
plot(nj.tree, tip.color = tip.color, cex = 0.7, main="NJ Tree Colored by Ploidy")
legend("topright", legend = levels(ploidy_factor), col = tip.cols, pch = 19)

############################################
# additional analyses using these distances
############################################
# 11
### create a dist object
colnames(lc.D.ind) <- rownames(lc.D.ind) 
lc.D.ind.dist<-as.dist(lc.D.ind, diag=T)
attr(lc.D.ind.dist, "Labels")<-rownames(lc.D.ind)          # name the rows of a matrix  

### plot NJ tree
par(cex = 0.4)
plot(nj(lc.D.ind))
write.tree(nj(lc.D.ind),file="NJ_lc_mixed_P.stampp.tree.tre")

#12
### heatmap of the pop distance matrix
pdf ("heatmap_Neis_Distances_lc_mixed_P.pdf", width=14, height=7)
heatmap.2(lc.D.ind, trace="none", cexRow=0.7, cexCol=0.7)
dev.off()

##########################################
### calculate AMOVA (differentiation among populations)
pops<-as.factor(substr(rownames(lc.D.ind),1,3))
(res <- pegas::amova(lc.D.ind.dist ~ pops))

### calculate AMOVA (differentiation among ploidy)
# extract ploidy for each individual
ploidy_factor <- as.factor(ploidy(lc.genlight))
table(ploidy_factor)

# if lc.D.ind is a dist object already, just use it
lc.D.ind.dist <- as.dist(lc.D.ind)

res_ploidy <- amova(lc.D.ind.dist ~ ploidy_factor)
res_ploidy


############
#   PCA 
############
#13
################# PCA based on real ploidy #######################
library(ggplot2)
library(dplyr)

# do PCA, retain first 300 axes (for later use in find.clusters)
pca.1 <- glPca(lc.genlight, nf=300) # see the modified function glPcaFast

# 1️⃣ Extract PCA coordinates
pca_df <- as.data.frame(pca.1$scores[, 1:3])
pca_df$Ploidy <- as.factor(ploidy(lc.genlight))

# 2️⃣ Ensure factor levels are explicitly set for all four ploidies
pca_df$Ploidy <- factor(pca_df$Ploidy, levels = c("2", "3", "4", "6"))  # update if you have other levels

# 3️⃣ Compute variance explained
pc1_var <- round((pca.1$eig[1] / sum(pca.1$eig)) * 100, 1)
pc2_var <- round((pca.1$eig[2] / sum(pca.1$eig)) * 100, 1)
pc3_var <- round((pca.1$eig[3] / sum(pca.1$eig)) * 100, 1)
pc4_var <- round((pca.1$eig[4] / sum(pca.1$eig)) * 100, 1)

# 4️⃣ Define colors for all ploidy levels
ploidy_colors <- c(
  "2" = "#8863ce",
  "3" = "#edf41b",
  "4" = "#33d5db",
  "6" = "#ec501b"
)

# 5️⃣ PCA plot
p_ploidy <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Ploidy)) +
  geom_point(size = 3, alpha = 0.85) +
  labs(
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC2 (", pc2_var, "%)"),
    title = "PCA of genotypes coloured by ploidy"
  ) +
  scale_color_manual(values = ploidy_colors, drop = FALSE) +  # drop=FALSE forces display of all
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )

# 6️⃣ Show and save
print(p_ploidy)
ggsave("PCA_by_ploidy_mix_ploidy_fixed_PC12.pdf", plot = p_ploidy, width = 8, height = 6)


############ PCA based on population locations #######################
#14
# 1️⃣ Define population from sample names (first 3 letters)
pop(lc.genlight) <- as.factor(substr(indNames(lc.genlight), 1, 3))
table(pop(lc.genlight))

# 2️⃣ Extract PCA coordinates
pca_df <- as.data.frame(pca.1$scores[, 1:3])
pca_df$Population <- factor(pop(lc.genlight))

# 3️⃣ Compute variance explained
pc1_var <- round((pca.1$eig[1] / sum(pca.1$eig)) * 100, 1)
pc2_var <- round((pca.1$eig[2] / sum(pca.1$eig)) * 100, 1)
pc3_var <- round((pca.1$eig[3] / sum(pca.1$eig)) * 100, 1)

# 4️⃣ Define colors for populations
pop_colors <- c(
  "GUD" = "#010102",  # red
  "BTR" = "#be4eb7",  # blue
  "JAI" = "#05720f",  # green
  "JOD" = "#aaf796",  # purple
  "SAT" = "#f1db14"   # orange
)

# 5️⃣ PCA plot
p_pop <- ggplot(pca_df, aes(x = PC2, y = PC3, color = Population)) +
  geom_point(size = 3, alpha = 0.85) +
  labs(
    x = paste0("PC2 (", pc2_var, "%)"),
    y = paste0("PC3 (", pc3_var, "%)"),
    title = "PCA of mixed-ploidy genotypes by population"
  ) +
  scale_color_manual(values = pop_colors, drop = FALSE) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )

# 6️⃣ Show and save
print(p_pop)
ggsave("PCA_by_population_mix_ploidy_PC23.pdf", plot = p_pop, width = 8, height = 6)



####################################
########### convert vcf to genodive 
####################################
#convert VCF file to get an input file for Genodive (up to 8-allels multiallelic)
#groups_file.txt is a file with two columns: fist with individual names, second with corresponding pop names
#usage: Rscript convert_genodive.R realdata.vcf.gz groups_file.txt output realdata_genodive.txt
#created by Joern Gerchen 2024-03-17
setwd("/Users/praveenp/Desktop/Chapter_2_analysis/27_10_25/Genodive analysis/")

library(vcfR)

args = commandArgs(trailingOnly=TRUE)
vcf_file<-args[1]
vcf_input<-read.vcfR("ddRAD_freebayse_mixed_ploidy_noIndel_biallelic_depth_GQ30_Q30_mac3_snp_maxM98.vcf")
population_file<-args[2]
output_file<-args[3]

#read population table
population_table<-read.table("sample_name_ploidy", header=FALSE, sep="\t", col.names = c("Individual", "Pop"))
n_populations=length(unique(population_table$Pop))
genotypes<-extract.gt(vcf_input, element="GT")

#replaces phase with regular genotypes
genotypes<-gsub("\\|", "/", genotypes)

#check table
#Encode genotypes numerically
allele_vector<-names((table(genotypes)))
allele_vector_sorted<-allele_vector[order(nchar(allele_vector), allele_vector, decreasing = TRUE)]
for(i in 1:length(allele_vector_sorted)){
  v_split<-unlist(strsplit(allele_vector_sorted[i], "/"))
  if(i==1){
    max_ploidy<-length(v_split)
  }
  v_inc<-sort(v_split)
  for(n in 1:length(v_split)){
    v_inc[n]<-as.character(as.numeric(v_split[n])+1)
  }
  i_inc<-paste(v_inc, collapse = "")
  genotypes<-gsub(allele_vector_sorted[i], i_inc, genotypes)
}

#Transpose genotypes
gen_transp<-t(genotypes)
sample_names<-rownames(gen_transp)
sample_pops<-sample_names

#Assign each individual to its population/ploidy
for(sample_name_i in 1:length(sample_names)){
  sample_pops[sample_name_i]<-as.character(population_table$Pop[population_table$Individual==sample_names[sample_name_i]])
}
loci_output<-c("pop", "ind", colnames(gen_transp))
sample_pops_numeric<-as.character(match(sample_pops, unique(population_table$Pop)))

#Build the Genodive output table
output_table<- cbind(as.matrix(sample_pops_numeric), as.matrix(sample_names), gen_transp)
n_individuals=ncol(genotypes)
n_loci<-nrow(genotypes)

#Write the Genodive file
#cat("Genodive input file generated from VCF\n", file=output_file)
cat("Genodive input file generated from VCF\n", file="genodive_out_all_pop_specific")
cat(paste(as.character(n_individuals),as.character(n_populations),as.character(n_loci),as.character(max_ploidy),"1\n", sep="\t"), file="genodive_out_all_pop_specific", append=TRUE)

for(pop in unique(population_table$Pop)){
  cat(paste(pop,"\n", sep=""), file="genodive_out_all_pop_specific", append=TRUE)
}

cat(paste(paste(loci_output, collapse = "\t"), "\n", sep=""), file="genodive_out_all_pop_specific", append=TRUE)
write.table(output_table, file="genodive_out_all_pop_specific", append=TRUE, quote=FALSE, na="0", sep="\t", row.names = FALSE, col.names = FALSE)




