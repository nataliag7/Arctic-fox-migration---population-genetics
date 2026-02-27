# 03_fst.R
# Pairwise FST between regions using SNPRelate::snpgdsFst

# BiocManager::install("SNPRelate")
# install.packages(c("tidyverse", "reshape2"))

library(SNPRelate)
library(gdsfmt)
library(tidyverse)
library(reshape2)

bed_file <- "data/AF.imputed.thin.bed"
bim_file <- "data/AF.imputed.thin.bim"
fam_file <- "data/AF.imputed.thin.fam"
gds_file <- "data/AF.imputed.thin.gds"

popinfo_file <- "data/sample_popinfo.tsv"

if (!dir.exists("figures")) dir.create("figures", showWarnings = FALSE)

#----------------------------------------------#
# Convert PLINK -> GDS if not already present  #
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

regions <- unique(popinfo$Region)

#-----------------------------------------------------------#
# Compute pairwise FST for each pair of regions             #
#-----------------------------------------------------------#

pairwise_fst_matrix <- function(group_vector, labels) {

  ng <- length(labels)
  mat <- matrix(NA_real_, nrow = ng, ncol = ng,
                dimnames = list(labels, labels))

  for (i in seq_len(ng)) {
    for (j in seq_len(ng)) {
      if (j <= i) next

      keep <- group_vector %in% c(labels[i], labels[j])

      fst_obj <- snpgdsFst(
        genofile,
        population   = factor(group_vector[keep]),
        sample.id    = sample_ids[keep],
        method       = "W&C84",
        autosome.only = TRUE,
        remove.monosnp = TRUE,
        verbose      = TRUE
      )

      # Overall Fst estimate
      mat[labels[i], labels[j]] <- fst_obj$Fst
      mat[labels[j], labels[i]] <- fst_obj$Fst
    }
  }

  diag(mat) <- 0
  mat
}

fst_regions <- pairwise_fst_matrix(popinfo$Region, regions)

#--------------------------------#
# Heatmap of FST (by regions)    #
#--------------------------------#

fst_df <- melt(fst_regions, varnames = c("Region1", "Region2"), value.name = "Fst")

p_fst <- ggplot(fst_df, aes(x = Region1, y = Region2, fill = Fst)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "magma", na.value = "grey95") +
  coord_equal() +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Pairwise FST between regions (Weir & Cockerham 1984)",
    x     = "",
    y     = "",
    fill  = "FST"
  )

ggsave("figures/fst_regions.png", p_fst, width = 5.5, height = 4.5, dpi = 300)

snpgdsClose(genofile)
