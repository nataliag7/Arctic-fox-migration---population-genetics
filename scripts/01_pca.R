# 01_pca.R
# PCA of Arctic fox SNP data using SNPRelate

# Install if needed:
# BiocManager::install("SNPRelate")
# install.packages(c("tidyverse"))

library(SNPRelate)
library(gdsfmt)
library(tidyverse)

#-------------#
# File paths  #
#-------------#

bed_file <- "data/AF.imputed.thin.bed"
bim_file <- "data/AF.imputed.thin.bim"
fam_file <- "data/AF.imputed.thin.fam"
gds_file <- "data/AF.imputed.thin.gds"

popinfo_file <- "data/sample_popinfo.tsv"

if (!dir.exists("figures")) dir.create("figures", showWarnings = FALSE)

#----------------------------------------------#
# Convert PLINK -> GDS (only first time)       #
#----------------------------------------------#

if (!file.exists(gds_file)) {
  snpgdsBED2GDS(
    bed.fn   = bed_file,
    bim.fn   = bim_file,
    fam.fn   = fam_file,
    out.gdsfn = gds_file,
    cvt.chr  = "char"
  )
}

genofile <- snpgdsOpen(gds_file)

#-------------------------#
# Sample info & metadata  #
#-------------------------#

sample_ids <- read.gdsn(index.gdsn(genofile, "sample.id"))

popinfo <- read_tsv(popinfo_file, col_types = cols())

# Keep only samples present in GDS and align order
popinfo <- popinfo %>%
  filter(Sample %in% sample_ids) %>%
  mutate(Sample = factor(Sample, levels = sample_ids)) %>%
  arrange(Sample)

#---------------------------#
# Run PCA with SNPRelate    #
#---------------------------#

pca <- snpgdsPCA(
  genofile,
  sample.id     = as.character(popinfo$Sample),
  autosome.only = TRUE,
  remove.monosnp = TRUE,
  maf           = 0.05,
  missing.rate  = 0.05,
  num.thread    = 2
)

pc_percent <- pca$varprop * 100

pca_df <- tibble(
  Sample = pca$sample.id,
  PC1    = pca$eigenvect[, 1],
  PC2    = pca$eigenvect[, 2],
  PC3    = pca$eigenvect[, 3]
) %>%
  left_join(popinfo, by = "Sample")

pc1_lab <- paste0("PC1 (", round(pc_percent[1], 1), "%)")
pc2_lab <- paste0("PC2 (", round(pc_percent[2], 1), "%)")
pc3_lab <- paste0("PC3 (", round(pc_percent[3], 1), "%)")

#---------------------------#
# Plots                     #
#---------------------------#

p_pc1_pc2 <- ggplot(pca_df, aes(x = PC1, y = PC2,
                                colour = Region, shape = Country)) +
  geom_point(size = 2, alpha = 0.9) +
  theme_bw(base_size = 12) +
  labs(
    title = "PCA – PC1 vs PC2",
    x = pc1_lab,
    y = pc2_lab
  )

p_pc2_pc3 <- ggplot(pca_df, aes(x = PC2, y = PC3,
                                colour = Region, shape = Country)) +
  geom_point(size = 2, alpha = 0.9) +
  theme_bw(base_size = 12) +
  labs(
    title = "PCA – PC2 vs PC3",
    x = pc2_lab,
    y = pc3_lab
  )

ggsave("figures/pca_pc1_pc2.png", p_pc1_pc2, width = 7, height = 5, dpi = 300)
ggsave("figures/pca_pc2_pc3.png", p_pc2_pc3, width = 7, height = 5, dpi = 300)

snpgdsClose(genofile)
