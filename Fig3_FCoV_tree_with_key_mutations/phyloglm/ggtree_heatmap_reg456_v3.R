library(ggtree)
library(ggplot2)
library(ape)
library(dplyr)

regions <- list(
  reg4 = list(
    tree = "FCoV1_S_all_seq_cleaned_cds_linsi_curated_fullllength_fix_dates_withSeCoV_reg4.fas.ingroup.treefile.nwk",
    data = "tree_annotation_biotype_SNP_reg4.csv"
  ),
  reg5 = list(
    tree = "FCoV1_S_all_seq_cleaned_cds_linsi_curated_fullllength_fix_dates_withSeCoV_reg5.fas.ingroup.treefile.nwk",
    data = "tree_annotation_biotype_SNP_reg5.csv"
  ),
  reg6 = list(
    tree = "FCoV1_S_all_seq_cleaned_cds_linsi_curated_fullllength_fix_dates_withSeCoV_reg6.fas.ingroup.treefile.nwk",
    data = "tree_annotation_biotype_SNP_reg6.csv"
  )
)


binarize_snp <- function(vec) {
  vec_clean <- vec
  vec_clean[vec_clean == "-"] <- NA
  tab <- sort(table(vec_clean), decreasing = TRUE)
  ref_allele <- names(tab)[1]
  ifelse(is.na(vec_clean), NA, ifelse(vec_clean == ref_allele, "WT", "Mut"))
}

build_tree_heatmap <- function(tree_file, data_file, region_label) {

  tr <- read.tree(tree_file)
  df <- read.csv(data_file, stringsAsFactors = FALSE)
  df <- df[match(tr$tip.label, df$Accession), ]
  rownames(df) <- df$Accession

  snp_cols <- grep("^SNP_", colnames(df), value = TRUE)

  ## heatmap matrix: rows = tip label, cols = SNP site, values "WT"/"Mut"
  mat <- as.data.frame(sapply(snp_cols, function(s) binarize_snp(df[[s]])))
  rownames(mat) <- rownames(df)
  colnames(mat) <- snp_cols

  ## tip metadata for biotype coloring, joined onto the tree with %<+%
  meta <- data.frame(label = tr$tip.label, biotype = df$biotype)

  p <- ggtree(tr, layout = "rectangular", linewidth = 0.5) %<+% meta +
    geom_tippoint(aes(color = biotype), size = 3) +
    scale_color_manual(values = c(FECV = "steelblue3", FIPV = "firebrick3"),
                        name = "Biotype") +
    theme(legend.position = "left")


  base_mapping <- p$mapping

  p2 <- gheatmap(
    p, mat,
    offset = max(node.depth.edgelength(tr)) * 0.05,
    width = 0.05 * ncol(mat),
    colnames_angle = 90, colnames_offset_y = 1, font.size = 2.5,
    color = "grey60"
  )


  attr(p2, "mapping") <- base_mapping

    p2 <- p2 +
    scale_fill_manual(values = c(WT = "grey85", Mut = "#0F007C"),
                       name = "Mutation state", na.value = "white") +
    ggtitle(paste0(region_label, ": tree (tip color = biotype) + candidate mutation heatmap")) +
    theme(plot.title = element_text(size = 12, face = "bold"))

  p2
}


for (region_name in names(regions)) {
  p <- build_tree_heatmap(
    regions[[region_name]]$tree,
    regions[[region_name]]$data,
    region_label = paste0("Region ", sub("reg", "", region_name))
  )
  ggsave(paste0("ggtree_heatmap_", region_name, ".pdf"), p, width = 10, height = 14)
  ggsave(paste0("ggtree_heatmap_", region_name, ".png"), p, width = 10, height = 14, dpi = 300)
}

## To combine all three regions into one multi-panel figure (e.g. for a
## manuscript Fig. 3), use patchwork:
##   library(patchwork)
##   p4 / p5 / p6   # or p4 | p5 | p6 for side-by-side
