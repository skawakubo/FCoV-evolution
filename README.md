# FCoV-evolution

Data and analysis code for:

Izumi Kida†, Shusuke Kawakubo†, Terza Brostoff, Petra Cerna, Michael R. Lappin,
Patricia A. Pesavento, Genotype to Phenotype Japan (G2P-Japan) Consortium,
Kei Sato, Keita Matsuno*.
**Genomic and functional insights into spike protein determinants driving
feline coronavirus biotype switching.**

† These authors contributed equally
\* Correspondence

## Repository structure

### `seqdata/` — sequence datasets used throughout the study
- `FCoV_wg_seqs/` — FCoV-1 whole-genome sequences, metadata, and the curation notebook
- `FCoV_spike_seqs/` — curated spike CDS/AA alignments, date-annotated metadata,
  CD-HIT-reduced sets (90/95/98%), and per-NRR region alignments (`NRR/`)
- `spike_domain0_del/` — alignments for the spike domain 0 deletion analysis

### `Fig1_FIPV_vs_FECV_Fishers_exact_test/` — identification of biotype-associated sites
Fisher's exact tests comparing FIPV and FECV sequences, per gene.
- `Spike/` — spike analyses (full length, furin cleavage site, M1058), volcano-plot
  inputs, per-residue composition plots, and Fig. 1D–F source data
- `ORFs/`, `3c/`, `7b/` — the same test applied to the remaining ORFs
- `extract_ORF_FCoV1.ipynb` — extracts each ORF from the whole-genome alignment

### `Fig2_FCoV_tree_with_key_mutations/` — phylogenetic analyses
- `spike_gard/` — GARD recombination screening and IQ-TREE trees for each
  non-recombinant region (NRR); `gard/` holds the whole-genome Alphacoronavirus-1 screen
- `treetime/` — TreeTime ancestral reconstruction per NRR and branch-mutation summaries
- `phyloglm/` — phylogenetic logistic regression of biotype against key mutations,
  plus the ggtree heatmap figures
- `mutation_combination_effect/` — co-occurrence and combined-effect analyses of key mutations
- `old_RDP_tree_file/` — earlier RDP-based analysis, retained for reference

### `Fig3_intrahost_mutations/` — intrahost (low-frequency) variants from NGS data
- `isnv_pileup.py`, `pileup_base_counts.py`, `base_composition.sh` — variant and
  base-composition calling from mapped BAMs
- `spike_overview.py` — summary figure across all samples along the spike CDS
- `refs/` — per-sample spike CDS references
- `plots_stats/` — per-sample outputs (regional counts, codon tables, plots)
- `examples/` — a worked example (`FECV3/`) including the variant-calling script
  and position bookkeeping

## Environment
`environment.yml` describes the conda environment used for the Python/Jupyter analyses:

```bash
conda env create -f environment.yml
conda activate FIP
