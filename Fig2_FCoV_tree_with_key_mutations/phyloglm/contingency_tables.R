regions <- list(
  reg4 = "tree_annotation_biotype_SNP_reg4.csv",
  reg5 = "tree_annotation_biotype_SNP_reg5.csv",
  reg6 = "tree_annotation_biotype_SNP_reg6.csv"
)

binarize_snp <- function(vec) {
  vec_clean <- vec
  vec_clean[vec_clean == "-"] <- NA
  tab <- sort(table(vec_clean), decreasing = TRUE)
  ref_allele <- names(tab)[1]
  bin <- ifelse(is.na(vec_clean), NA, ifelse(vec_clean == ref_allele, 0, 1))
  list(bin = bin, ref = ref_allele, alt = names(tab)[-1])
}

long_rows <- list()
txt_out <- character()

for (region_name in names(regions)) {
  df <- read.csv(regions[[region_name]], stringsAsFactors = FALSE)
  snp_cols <- grep("^SNP_", colnames(df), value = TRUE)

  for (site in snp_cols) {
    b <- binarize_snp(df[[site]])
    mut_state <- factor(b$bin, levels = c(0, 1),
                         labels = c(paste0(b$ref, " (ref)"),
                                    paste0(paste(b$alt, collapse = "/"), " (alt)")))
    tab2x2 <- table(biotype = df$biotype, mutation = mut_state, useNA = "no")

    txt_out <- c(txt_out,
                 sprintf("== %s / %s (ref=%s, alt=%s) ==", region_name, site, b$ref,
                         paste(b$alt, collapse = "/")),
                 capture.output(print(tab2x2)),
                 "")

    df_tab <- as.data.frame(tab2x2)
    df_tab$region <- region_name
    df_tab$site <- site
    long_rows[[paste(region_name, site)]] <- df_tab[, c("region", "site", "biotype", "mutation", "Freq")]
  }
}

long_table <- do.call(rbind, long_rows)
rownames(long_table) <- NULL

write.csv(long_table, "/contingency_tables_reg4_5_6_long.csv", row.names = FALSE)
writeLines(txt_out, "ontingency_tables_reg4_5_6.txt")
