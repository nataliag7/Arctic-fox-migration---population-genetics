# 04_nj_tree.R
# Neighbor-joining tree based on genome-wide IBS

# BiocManager::install("SNPRelate")
# install.packages(c("tidyverse", "ape", "RColorBrewer"))

library(SNPRelate)
library(gdsfmt)
library(tidyverse)
library(ape)
library(RColorBrewer)

bed_file <- "data/AF.imputed.thin.bed"
bim_file <- "data/AF.imputed.thin.bim"
fam_file <- "data/AF.imputed.thin.fam"
gds_file <- "data/AF.imputed.thin.gds"

popinfo_file <- "data/sample_popinfo.tsv"

if (!dir.exists("figures")) dir.create("figures", showWarnings = FALSE)

#----------------------------------------------#
# Convert PLINK -> GDS if needed               #
#----------------------------------------------#

if (!file.exists(gds_file)) {
  snpgdsBED2GDS(
    bed.fn    = bed_file,
    bim.fn    = bim_file,
    fam.fn    = fam_file,
    out.gdsfn = gds_file,
    cvt.chr   = "char"
  )
}

genofile   <- snpgdsOpen(gds_file)
sample_ids <- read.gdsn(index.gdsn(genofile, "sample.id"))

popinfo <- read_tsv(popinfo_file, col_types = cols()) %>%
  filter(Sample %in% sample_ids) %>%
  mutate(Sample = factor(Sample, levels = sample_ids)) %>%
  arrange(Sample)

#----------------------------#
# Identity-by-state (IBS)    #
#----------------------------#

ibs <- snpgdsIBS(
  genofile,
  sample.id     = as.character(popinfo$Sample),
  autosome.only = TRUE,
  remove.monosnp = TRUE,
  maf           = 0.05,
  missing.rate  = 0.05,
  num.thread    = 2,
  verbose       = TRUE
)

# Convert IBS proportion (0–1) to distance
dist_mat <- 1 - ibs$ibs
rownames(dist_mat) <- colnames(dist_mat) <- as.character(popinfo$Sample)

#----------------------------#
# Neighbor-joining tree      #
#----------------------------#

nj_tree <- nj(as.dist(dist_mat))

# Map tip labels back to metadata order
tip_regions <- popinfo$Region[match(nj_tree$tip.label, popinfo$Sample)]

region_levels <- unique(popinfo$Region)
region_cols   <- setNames(
  brewer.pal(max(3, min(8, length(region_levels))), "Set1"),
  region_levels
)

tip_cols <- region_cols[tip_regions]

png("figures/nj_tree.png", width = 800, height = 800)
plot(
  nj_tree,
  tip.color = tip_cols,
  cex       = 0.8,
  main      = "Neighbor-joining tree (IBS distance)"
)
legend(
  "topleft",
  legend = names(region_cols),
  col    = region_cols,
  pch    = 19,
  cex    = 0.8,
  bty    = "n"
)
dev.off()

snpgdsClose(genofile)
