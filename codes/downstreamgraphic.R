library(ggplot2)
library(dplyr)

# Read predicted score data and standardize column names
solo_scores <- read.csv("E:/CMML/part2/ICA/paper_result/solo/solo_doublet_scores.csv") %>%
  dplyr::rename(predicted_score = predicted_doublet) %>%   # Rename column using dplyr's rename
  dplyr::mutate(Method = "Solo")

doublet_scores <- read.csv("E:/CMML/part2/ICA/paper_result/doubletfinder/doublet_parameters.csv") %>%
  mutate(Method = "DoubletFinder")  # Already has predicted_score

scrublet_scores <- read.csv("E:/CMML/part2/ICA/paper_result/scrublet/scrublet_results.csv",
                            header = TRUE, stringsAsFactors = FALSE) %>%
  dplyr::rename(predicted_score = scrublet_score) %>%
  dplyr::mutate(Method = "Scrublet")

# Combine data
scores_df <- bind_rows(
  solo_scores %>% select(predicted_score, Method),
  doublet_scores %>% select(predicted_score, Method),
  scrublet_scores %>% select(predicted_score, Method)
)

# Plot the density curves comparing the predicted score distributions from three methods
p_density <- ggplot(scores_df, aes(x = predicted_score, fill = Method, color = Method)) +
  geom_density(alpha = 0.3, bw = 0.005) +      # Adjust bandwidth for more detail
  scale_x_continuous(limits = c(0, 0.25)) +      # Set x-axis range
  labs(title = "Comparison of Predicted Doublet Scores Density Distribution", 
       x = "Predicted Score", 
       y = "Density") +
  theme_minimal()
ggsave("E:/CMML/part2/ICA/paper_result/combined_predicted_score_density.png", 
       plot = p_density, dpi = 300)


library(ggplot2)
library(dplyr)

# Read doublet proportion data for each method by cluster
solo_cluster <- read.csv("E:/CMML/part2/ICA/paper_result/solo/solo_cluster_doublet_proportions.csv") %>%
  mutate(Method = "Solo")
doublet_cluster <- read.csv("E:/CMML/part2/ICA/paper_result/doubletfinder/doubletfinder_cluster_doublet_proportions.csv") %>%
  mutate(Method = "DoubletFinder")
scrublet_cluster <- read.csv("E:/CMML/part2/ICA/paper_result/scrublet/scrublet_cluster_doublet_proportions.csv") %>%
  mutate(Method = "Scrublet")

# Combine cluster data
cluster_df <- bind_rows(solo_cluster, doublet_cluster, scrublet_cluster)

# Plot bar chart showing the doublet rate by cluster and method
p_cluster <- ggplot(cluster_df, aes(x = as.factor(seurat_clusters), y = Doublet_Rate, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Comparison of Doublet Proportions Across Clusters", 
       x = "Cluster", 
       y = "Doublet Rate") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggsave("E:/CMML/part2/ICA/paper_result/combined_cluster_doublet_proportions.png", 
       plot = p_cluster, dpi = 300)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(grid)

# Read DE results
solo_de <- read.csv("E:/CMML/part2/ICA/paper_result/solo/solo_DE_results.csv", row.names = 1)
doublet_de <- read.csv("E:/CMML/part2/ICA/paper_result/doubletfinder/doubletfinder_DE_results.csv", row.names = 1)
scrublet_de <- read.csv("E:/CMML/part2/ICA/paper_result/scrublet/scrublet_DE_results.csv", row.names = 1)

# Calculate -log10(p_val_adj) and save row names as 'gene'
solo_de$negLog10_p <- -log10(solo_de$p_val_adj)
solo_de$gene <- rownames(solo_de)

doublet_de$negLog10_p <- -log10(doublet_de$p_val_adj)
doublet_de$gene <- rownames(doublet_de)

scrublet_de$negLog10_p <- -log10(scrublet_de$p_val_adj)
scrublet_de$gene <- rownames(scrublet_de)

# Set thresholds (use abs(log2 Fold Change) > 1 and -log10(p_val_adj) > 2 as examples)
fold_change_threshold <- 1
pvalue_threshold <- 2   # corresponding to p_val_adj < 0.01

# Volcano plot for Solo (set points to blue)
p1 <- ggplot(solo_de, aes(x = avg_log2FC, y = negLog10_p)) +
  geom_point(alpha = 0.6, color = "blue") +
  geom_vline(xintercept = c(-fold_change_threshold, fold_change_threshold),
             linetype = "dashed", color = "grey") +
  geom_hline(yintercept = pvalue_threshold,
             linetype = "dashed", color = "grey") +
  geom_text_repel(data = subset(solo_de,
                                negLog10_p > pvalue_threshold & abs(avg_log2FC) > fold_change_threshold),
                  aes(label = gene),
                  size = 3,
                  max.overlaps = 10,
                  box.padding = 0.5,
                  point.padding = 0.5,
                  segment.size = 0.3) +
  theme_minimal() +
  labs(title = "Solo DE Volcano",
       x = "log2 Fold Change",
       y = "-log10(Adjusted p-value)")

# Volcano plot for DoubletFinder (set points to red)
p2 <- ggplot(doublet_de, aes(x = avg_log2FC, y = negLog10_p)) +
  geom_point(alpha = 0.6, color = "red") +
  geom_vline(xintercept = c(-fold_change_threshold, fold_change_threshold),
             linetype = "dashed", color = "grey") +
  geom_hline(yintercept = pvalue_threshold,
             linetype = "dashed", color = "grey") +
  geom_text_repel(
    data = subset(doublet_de, negLog10_p > pvalue_threshold & abs(avg_log2FC) > fold_change_threshold),
    aes(label = gene),
    size = 3,
    max.overlaps = 10,
    box.padding = 0.5,
    point.padding = 0.5,
    segment.size = 0.3
  ) +
  theme_minimal() +
  labs(title = "DoubletFinder DE Volcano",
       x = "log2 Fold Change",
       y = "-log10(Adjusted p-value)")

# Volcano plot for Scrublet (set points to green)
p3 <- ggplot(scrublet_de, aes(x = avg_log2FC, y = negLog10_p)) +
  geom_point(alpha = 0.6, color = "green") +
  geom_vline(xintercept = c(-fold_change_threshold, fold_change_threshold),
             linetype = "dashed", color = "grey") +
  geom_hline(yintercept = pvalue_threshold,
             linetype = "dashed", color = "grey") +
  geom_text_repel(data = subset(scrublet_de,
                                negLog10_p > pvalue_threshold & abs(avg_log2FC) > fold_change_threshold),
                  aes(label = gene),
                  size = 3,
                  max.overlaps = 10,
                  box.padding = 0.5,
                  point.padding = 0.5,
                  segment.size = 0.3) +
  theme_minimal() +
  labs(title = "Scrublet DE Volcano",
       x = "log2 Fold Change",
       y = "-log10(Adjusted p-value)")

# Arrange the three volcano plots side by side

combined_grob <- arrangeGrob(p1, p2, p3, ncol = 3)
# Save the combined volcano plots figure
ggsave("E:/CMML/part2/ICA/paper_result/combined_DE_volcano_plots.png",
       plot = combined_grob,
       width = 15,   # 根据需要调整宽度
       height = 5,   # 根据需要调整高度
       dpi = 300)
library(ggplot2)
library(dplyr)

# 读取 CSV 文件，并重命名 UMAP 坐标列（统一为 UMAP1 和 UMAP2）
doublet_umap <- read.csv("E:/CMML/part2/ICA/paper_result/doubletfinder/doublet_umap.csv") %>%
      rename(UMAP1 = umap_1, UMAP2 = umap_2)

# 绘制 UMAP 散点图，采用固定颜色
p_umap <- ggplot(doublet_umap, aes(x = UMAP1, y = UMAP2)) +
  geom_point(alpha = 0.7, size = 1.5, color = "blue") +
  labs(title = "DoubletFinder UMAP Plot") +
  theme_minimal()

# 显示图形
print(p_umap)

# 保存图形
ggsave("E:/CMML/part2/ICA/paper_result/doubletfinder/umap_plot.png", 
       plot = p_umap, width = 8, height = 6, dpi = 300)