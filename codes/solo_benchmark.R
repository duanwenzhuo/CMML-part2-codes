# ========= Load required packages =========
library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(slingshot)
library(data.table)

# ========= 1. Read CSV Files =========
# Read metadata and prediction scores (using SOLO version as an example)
solo_obs <- read.csv("E:/CMML/part2/ICA/paper_result/solo/solo_obs_results.csv", row.names = 1)
solo_perf <- read.csv("E:/CMML/part2/ICA/paper_result/solo/solo_doublet_scores.csv")
print("SOLO performance metrics:")
print(solo_perf)

# ========= 2. Create the Seurat Object =========
# Read raw count data (rows: genes, columns: cells)
count_matrix <- readRDS("E:/CMML/part2/ICA/paper_result/counts_matrix.rds")
if (is.null(rownames(count_matrix))) {
  rownames(count_matrix) <- paste0("gene_", 1:nrow(count_matrix))
}
if (is.null(colnames(count_matrix))) {
  colnames(count_matrix) <- paste0("cell_", 1:ncol(count_matrix))
}
# Create the Seurat object using the count matrix and metadata
seurat_obj <- CreateSeuratObject(counts = count_matrix, meta.data = solo_obs)

# Data preprocessing: normalization, variable features selection, scaling, and PCA
seurat_obj <- NormalizeData(seurat_obj) %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA(seurat_obj, verbose = FALSE)

# Set the 'cell_annotation' from metadata as the true labels (stored in new column 'doublet_label')
seurat_obj@meta.data$doublet_label <- as.factor(solo_obs$cell_annotation)
Idents(seurat_obj) <- seurat_obj@meta.data$doublet_label

# ========= 3. Differential Expression (DE) Analysis =========
# Assuming a comparison between "doublet" and "singlet"
de_results <- FindMarkers(seurat_obj, ident.1 = "doublet", ident.2 = "singlet")
write.csv(de_results, "E:/CMML/part2/ICA/paper_result/solo/solo_DE_results.csv", row.names = TRUE)
print("Differential expression analysis completed and saved to solo_DE_results.csv")

# ========= 4. Clustering and UMAP Visualization =========
# Perform clustering on cells and compute UMAP (using PCA results, dims 1:10)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)

# Plot UMAP showing clusters
DimPlot(seurat_obj, reduction = "umap", label = TRUE) + ggtitle("Cluster Results")

# Calculate the proportion of true labels (doublet_label) in each cluster
cluster_doublet_prop <- seurat_obj@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(Doublet_Rate = mean(as.numeric(doublet_label)))
print(cluster_doublet_prop)
write.csv(cluster_doublet_prop, 
          "E:/CMML/part2/ICA/paper_result/solo/solo_cluster_doublet_proportions.csv", 
          row.names = FALSE)

# ========= 5. Trajectory Analysis =========
# Based on Seurat's UMAP results and clustering, use Slingshot to infer cell trajectories
umap_coords <- Embeddings(seurat_obj, "umap")
cluster_labels <- seurat_obj@meta.data$seurat_clusters
doublet_label <- seurat_obj@meta.data$doublet_label
# Color cells based on true labels: mark "doublet" in red, others in blue
cell_colors <- ifelse(doublet_label == "doublet", "red", "blue")
# Select the cluster with the most cells as the starting cluster
init_clust <- names(which.max(table(cluster_labels)))
sling_obj <- slingshot(umap_coords, clusterLabels = cluster_labels, start.clus = init_clust)

# Plot UMAP with trajectory overlay
plot(umap_coords, col = cell_colors, pch = 16, asp = 1,
     main = "Trajectory Analysis: Distribution of Doublet Cells (Slingshot)")
curves_list <- slingCurves(sling_obj)
for (i in seq_along(curves_list)) {
  curve_i <- curves_list[[i]]$curve
  if (!is.null(curve_i)) {
    lines(curve_i[, 1], curve_i[, 2], col = "black", lwd = 2)
  }
}
legend("topright", legend = c("Doublet", "Singlet"), col = c("red", "blue"), pch = 16)

# ========= 6. UMAP Visualization Based on Predicted Scores =========
# Use FeaturePlot to display the predicted doublet score ("predicted_doublet") from metadata
FeaturePlot(seurat_obj, features = "predicted_doublet") +
  ggtitle("UMAP - Predicted Doublet Scores")
#得到corrected data 
# 假设两个文件中的细胞 ID 一致（这里利用行名或一个公共列进行合并）
solo_merged <- merge(solo_obs, solo_scores, by = "row.names", all.x = TRUE)
rownames(solo_merged) <- solo_merged$Row.names
solo_merged <- solo_merged[,-1]  # 移除多余的第一列

# 查看预测得分分布情况
summary(solo_merged$predicted_score)

# 设定阈值：例如去除得分最高的 40% 细胞，则保留得分低于 60% 分位数的细胞
threshold <- quantile(solo_merged$predicted_score, 0.6, na.rm = TRUE)
cat("使用的阈值为：", threshold, "\n")

# 筛选出预测得分低于该阈值的细胞，即认为为单细胞的 Corrected data
corrected_data <- solo_merged[solo_merged$predicted_score <= threshold, ]
cat("Corrected data 细胞数：", nrow(corrected_data), "\n")

# 如有需要，可导出 Corrected data 供后续差异表达分析等使用
write.csv(corrected_data, "E:/CMML/part2/ICA/paper_result/solo/solo_corrected_data.csv", row.names = TRUE)
