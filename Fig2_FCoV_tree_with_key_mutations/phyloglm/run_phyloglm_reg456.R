library(ape)
library(phylolm)

set.seed(1)


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

## Optional: specify a known causal allele per SNP site to binarize
## "that allele" vs "everything else" instead of "majority vs pooled minor".
## Leave empty to use the default (majority vs pooled minor) for all sites.
## Example: known_causal_allele <- c(SNP_814 = "S")
known_causal_allele <- c()

## ---- 2. Helper: binarize one multi-allelic SNP column -----------------------

binarize_snp <- function(vec, causal_allele = NULL) {
  vec_clean <- vec
  vec_clean[vec_clean == "-"] <- NA

  if (!is.null(causal_allele)) {
    bin <- ifelse(is.na(vec_clean), NA, ifelse(vec_clean == causal_allele, 1, 0))
    tab <- table(vec_clean)
    return(list(bin = bin, ref = paste0("not-", causal_allele),
                alt = causal_allele, alt_counts = tab[causal_allele]))
  }

  tab <- sort(table(vec_clean), decreasing = TRUE)
  ref_allele <- names(tab)[1]
  bin <- ifelse(is.na(vec_clean), NA, ifelse(vec_clean == ref_allele, 0, 1))
  list(bin = bin, ref = ref_allele, alt = names(tab)[-1], alt_counts = tab[-1])
}

## ---- 3. Helper: run phyloglm for one site within one region ------------------

run_site_phyloglm <- function(site, df, tr) {
  causal <- known_causal_allele[site]
  b <- binarize_snp(df[[site]], causal_allele = if (length(causal) == 1 && !is.na(causal)) causal else NULL)

  keep <- !is.na(b$bin)
  tr_i <- tr
  dropped_n <- sum(!keep)
  if (dropped_n > 0) tr_i <- drop.tip(tr, rownames(df)[!keep])

  mut01 <- setNames(b$bin[keep], rownames(df)[keep])[tr_i$tip.label]
  bio01 <- setNames(ifelse(df$biotype == "FIPV", 1, 0), rownames(df))[tr_i$tip.label]

  n_used  <- length(mut01)
  n_minor <- sum(mut01 == 1)

  base_row <- data.frame(
    site = site, ref_allele = b$ref,
    alt_alleles = paste(names(b$alt_counts), collapse = "/"),
    n_used = n_used, n_minor = n_minor,
    phyloglm_beta = NA, phyloglm_se = NA, phyloglm_OR = NA,
    phyloglm_p = NA, note = ""
  )

  if (n_minor < 2 || n_minor > (n_used - 2)) {
    base_row$note <- "skipped: minor allele count too extreme (<2 or >n-2)"
    return(base_row)
  }

  dat <- data.frame(biotype01 = bio01, mutation01 = as.numeric(mut01),
                     row.names = names(bio01))

  fit <- tryCatch(
    phyloglm(biotype01 ~ mutation01, data = dat, phy = tr_i,
             method = "logistic_MPLE", boot = 200),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    base_row$note <- paste("phyloglm error:", conditionMessage(fit))
    return(base_row)
  }

  s <- summary(fit)
  base_row$phyloglm_beta <- s$coefficients["mutation01", "Estimate"]
  base_row$phyloglm_se   <- s$coefficients["mutation01", "StdErr"]
  base_row$phyloglm_OR   <- exp(base_row$phyloglm_beta)
  base_row$phyloglm_p    <- s$coefficients["mutation01", "p.value"]
  base_row
}

## ---- 4. Loop over regions and sites -------------------------------------------

all_results <- list()

for (region_name in names(regions)) {
  cat("\n############ Region:", region_name, "############\n")

  tr <- read.tree(regions[[region_name]]$tree)
  df <- read.csv(regions[[region_name]]$data, stringsAsFactors = FALSE)
  stopifnot(setequal(tr$tip.label, df$Accession))
  df <- df[match(tr$tip.label, df$Accession), ]
  rownames(df) <- df$Accession
  tr$edge.length[tr$edge.length <= 0] <- 1e-7

  cat("  n tips:", Ntip(tr), " | rooted:", is.rooted(tr),
      " | binary:", is.binary(tr), " | zero-edges nudged:",
      sum(tr$edge.length == 1e-7), "\n")

  snp_cols <- grep("^SNP_", colnames(df), value = TRUE)

  for (site in snp_cols) {
    row <- run_site_phyloglm(site, df, tr)
    row <- cbind(region = region_name, row)
    cat(sprintf("  %-10s ref=%-3s alt=%-15s n=%3d n_minor=%3d  OR=%7s  p=%8s  %s\n",
                site, row$ref_allele, row$alt_alleles, row$n_used, row$n_minor,
                ifelse(is.na(row$phyloglm_OR), "NA", signif(row$phyloglm_OR, 3)),
                ifelse(is.na(row$phyloglm_p), "NA", signif(row$phyloglm_p, 3)),
                row$note))
    all_results[[paste(region_name, site)]] <- row
  }
}

## ---- 5. Combined summary table + multiple-testing correction -----------------

summary_table <- do.call(rbind, all_results)
rownames(summary_table) <- NULL

n_tests <- sum(!is.na(summary_table$phyloglm_p))
summary_table$phyloglm_p_bonferroni <- pmin(summary_table$phyloglm_p * n_tests, 1)
summary_table$phyloglm_p_BH <- p.adjust(summary_table$phyloglm_p, method = "BH")

summary_table <- summary_table[order(summary_table$phyloglm_p), ]

cat("\n\n===== COMBINED SUMMARY: reg4 + reg5 + reg6 (", n_tests, "tests) =====\n")
print(summary_table, digits = 3, row.names = FALSE)

out_path <- "phyloglm_reg4_5_6_results.csv"
write.csv(summary_table, out_path, row.names = FALSE)
cat("\nSaved:", out_path, "\n")
