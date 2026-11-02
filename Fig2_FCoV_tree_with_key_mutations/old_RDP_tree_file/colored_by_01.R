# ----------------------------------------------------------------              
# 1. 解析対象の設定（ここだけ書き換えればOK！）
# ----------------------------------------------------------------
# 解析したいポジションを指定
target_pos <- 811  

# ヒートマップ表示用の列名を生成（例：SNP_811）
target_col <- paste0("SNP_", target_pos)

# ----------------------------------------------------------------              
# 2. ライブラリの読み込み
# ----------------------------------------------------------------
library(ggplot2)
library(ggtree)
library(ggnewscale)
library(ape)
library(dplyr)
library(scales)

# ----------------------------------------------------------------              
# 3. データの読み込み
# ----------------------------------------------------------------
tree <- read.nexus('2344_annotated_tree.nexus')

# メタデータの読み込み
file_path <- '/Users/izumikida/FCoV-evolution/FCoV_tree_with_key_mutations/tree_annotation_biotype_SNP.csv'
temp_data <- read.csv(file_path, nrows = 1, check.names = FALSE)
col_types <- sapply(temp_data, function(x) "character") 
meta_data <- read.csv(file_path, colClasses = col_types, check.names = FALSE)

# Mutationデータの読み込み
mut_file_path <- '/Users/izumikida/FCoV-evolution/FCoV_tree_with_key_mutations/combine_NRR_branch_mutation_summary.tsv'
mutation_data <- read.delim(mut_file_path, colClasses = "character", check.names = FALSE)
mutation_data$pos <- as.numeric(mutation_data$pos)

# ----------------------------------------------------------------              
# 4. 系統樹の基本描図（FIPV=水色, FECV=灰色）
# ----------------------------------------------------------------
p <- ggtree(tree, layout = "rectangular") %<+% meta_data

p_final <- p + 
  geom_tippoint(aes(color = biotype), size = 4, alpha = 0.5) + 
  scale_color_manual(values = c("FECV" = "gray40", "FIPV" = "#8CA8FF")) + 
  geom_treescale(width = 0.1, x = 0, y = -5, fontsize = 2, offset = 1) +
  theme(legend.position = "right")

# --- (前略：1〜5までは変更なし) ---

# ----------------------------------------------------------------              
# 5. 変異獲得部位のマーク（Internal Node / Tip）
# ----------------------------------------------------------------
target_nodes <- mutation_data %>%
  filter(pos == target_pos) %>%
  pull(node) %>%
  unique()

target_df <- p$data %>%
  filter(label %in% target_nodes | node %in% target_nodes)

p_mutation_marked <- p_final + 
  geom_point(data = target_df, 
             aes(x = x, y = y, fill = isTip, shape = isTip), 
             size = 2.5, color = "gray10", stroke = 0.4, alpha = 0.7) + 
  scale_fill_manual(values = c("TRUE" = "darkorange", "FALSE" = "#006400"),   
                    labels = c("TRUE" = "External (Tip)", "FALSE" = "Internal Node"),
                    name = "Mutation Location") +
  scale_shape_manual(values = c("TRUE" = 21, "FALSE" = 24),
                     labels = c("TRUE" = "External (Tip)", "FALSE" = "Internal Node"),
                     name = "Mutation Location") +
  labs(title = paste("Acquisition of Mutation at Position", target_pos))

# ----------------------------------------------------------------              
# 6. ヒートマップの追加 (R = Gray, Mutations = Dark Purple)
# ----------------------------------------------------------------
# fillスケールをリセット
p_for_heatmap <- p_mutation_marked + new_scale_fill() 

data_heatmap <- meta_data[, target_col, drop = FALSE]
rownames(data_heatmap) <- meta_data$Accession

# データ内に存在する全アミノ酸を取得して色を割り当てる
present_aminos <- unique(unlist(data_heatmap))

# カラーパレットの作成
# 基本をすべて「濃い紫 (#4B0082)」にし、R だけを「灰色 (gray80)」にする
custom_colors <- setNames(rep("#4B0082", length(present_aminos)), present_aminos)

if("R" %in% names(custom_colors)) {
  custom_colors["R"] <- "gray80" # コンセンサスRを明るい灰色に
}

# ヒートマップの描画
final_plot <- gheatmap(p = p_for_heatmap, 
                       data = data_heatmap, 
                       offset = 0.08,       # Tip Pointとの距離を調整
                       width = 0.05, 
                       color = "white",     # タイルの枠線
                       colnames = TRUE, 
                       colnames_position = "top", 
                       font.size = 3) +
  # 作成したカスタムパレットを適用
  scale_fill_manual(values = custom_colors, 
                    name = paste0("Residue at ", target_pos), 
                    na.value = "white") +
  vexpand(0.1)

# ----------------------------------------------------------------              
# 7. 表示
# ----------------------------------------------------------------
print(final_plot)