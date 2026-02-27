# 02_admixture_admixplots.R
# Plot ADMIXTURE Q-matrices for different K

# install.packages(c("tidyverse", "RColorBrewer"))

library(tidyverse)
library(RColorBrewer)

popinfo_file <- "data/sample_popinfo.tsv"
admix_dir    <- "data/admixture"

if (!dir.exists("figures")) dir.create("figures", showWarnings = FALSE)

#------------------------#
# Sample metadata        #
#------------------------#

popinfo <- read_tsv(popinfo_file, col_types = cols())

# Order regions as they should appear in barplots
popinfo <- popinfo %>%
  mutate(
    Region  = factor(Region),
    Country = factor(Country)
  )

#--------------------------------------------------#
# Helper to plot ADMIXTURE for a single K          #
#--------------------------------------------------#

plot_admixture_for_K <- function(K,
                                 q_file  = file.path(admix_dir, sprintf("AF_K%d_best.Q", K)),
                                 out_png = file.path("figures", sprintf("admixture_k%d.png", K))) {
  if (!file.exists(q_file)) {
    stop("Q file not found: ", q_file)
  }

  qmat <- as.matrix(read_table(q_file, col_names = FALSE))
  colnames(qmat) <- paste0("Cluster_", seq_len(ncol(qmat)))

  if (nrow(qmat) != nrow(popinfo)) {
    stop("Number of rows in Q matrix (", nrow(qmat),
         ") does not match nrow(popinfo) (", nrow(popinfo), "). ",
         "Make sure popinfo is in the same individual order as the PLINK/ADMIXTURE input.")
  }

  df <- bind_cols(popinfo, as_tibble(qmat)) %>%
    mutate(
      # keep current order -> x axis will group by region via facet
      Sample = factor(Sample, levels = Sample)
    ) %>%
    pivot_longer(
      starts_with("Cluster_"),
      names_to  = "Cluster",
      values_to = "Ancestry"
    )

  nK <- length(unique(df$Cluster))
  palette_cols <- brewer.pal(max(3, min(8, nK)), "Set3")

  p <- ggplot(df, aes(x = Sample, y = Ancestry, fill = Cluster)) +
    geom_bar(stat = "identity", width = 1) +
    facet_grid(~ Region, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = palette_cols) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      panel.spacing = unit(0.2, "lines"),
      strip.background = element_rect(fill = "grey90"),
      strip.text.x = element_text(size = 9)
    ) +
    labs(
      title = paste("ADMIXTURE – K =", K),
      x     = "Individuals (grouped by region)",
      y     = "Ancestry proportion"
    )

  ggsave(out_png, p, width = 9, height = 4, dpi = 300)
  message("Saved: ", out_png)
}

#----------------------#
# Run for chosen Ks    #
#----------------------#

K_values <- c(5, 6, 7)

for (K in K_values) {
  plot_admixture_for_K(K)
}
