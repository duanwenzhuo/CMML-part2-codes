# ========= Load Required Packages =========
library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(slingshot)
library(PRROC)
library(data.table)

# ========= 1. Read CSV Files – Predicted Scores, Metadata, and Performance =========
# Using the SCRUBLET version (change file prefixes accordingly)
scrublet_obs <- read.csv("E:/CMML/part2/ICA/paper_result/scrublet/scrublet_obs_results.csv")
scrublet_perf <- read.csv("E:/CMML/part2/ICA/paper_result/scrublet/scrublet_results.csv")
# Optionally, you can also read threshold analysis if available:
# scrublet_thresh <- read.csv("E:/CMML/part2/ICA/paper_result/doubletfinder/doublet_threshold_analysis.csv")

print("SCRUBLET performance metrics:")
print(scrublet_perf)

# ========= 2. Differential Expression (DE) Analysis =========
# For comparing doublets and singlets based on true labels (assumed stored in "doublet_label")
# and predicted scores in "predicted_doublet".
# Here we load the raw count data (rows: genes; columns: cells) from an RDS file.
count_matrix <- readRDS("E:/CMML/part2/ICA/paper_result/counts_matrix.rds")
if (is.null(rownames(count_matrix))) {
  rownames(count_matrix) <- paste0("gene_", 1:nrow(count_matrix))
}
if (is.null(colnames(count_matrix))) {
  colnames(count_matrix) <- paste0("cell_", 1:ncol(count_matrix))
}
# Create the Seurat object by combining the count matrix with metadata.
seurat_obj <- CreateSeuratObject(counts = count_matrix, meta.data = scrublet_obs)

# Preprocess the data: normalization, variable gene selection, scaling, and PCA.
seurat_obj <- NormalizeData(seurat_obj) %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA(seurat_obj, verbose = FALSE)

# Use the metadata "cell_annotation" as the true label, stored in a new column "doublet_label"
labels <- scrublet_obs$cell_annotation
labels <- as.factor(labels)
seurat_obj@meta.data$doublet_label <- labels
Idents(seurat_obj) <- seurat_obj@meta.data$doublet_label

# Find genes differentially expressed between "doublet" and "singlet" cells.
de_results <- FindMarkers(seurat_obj, ident.1 = "doublet", ident.2 = "singlet")
write.csv(de_results, "E:/CMML/part2/ICA/paper_result/scrublet/scrublet_DE_results.csv", row.names = TRUE)
print("Differential expression analysis completed and saved to scrublet_DE_results.csv")

# ========= 3. Clustering and UMAP Visualization =========
# Perform clustering and compute UMAP using the PCA results (dims 1:10).
seurat_obj <- seurat_obj %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(verbose = FALSE)


seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)

# Plot the UMAP with cluster labels.
DimPlot(seurat_obj, reduction = "umap", label = TRUE) + ggtitle("Cluster Results")

# Calculate the proportion of doublets (from the true labels) in each cluster.
cluster_doublet_prop <- seurat_obj@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(Doublet_Rate = mean(as.numeric(doublet_label)))
print(cluster_doublet_prop)
write.csv(cluster_doublet_prop, 
          "E:/CMML/part2/ICA/paper_result/scrublet/scrublet_cluster_doublet_proportions.csv", 
          row.names = FALSE)

# ========= 4. Distributed Computing & Scalability Testing =========
# If Python generated a runtime and memory usage summary (runtime_memory_summary.csv),
# load the file to plot and compare performance.
if(file.exists("E:/CMML/part2/ICA/paper_result/scrublet/runtime_memory_summary.csv")){
  runtime_data <- fread("E:/CMML/part2/ICA/paper_result/scrublet/runtime_memory_summary.csv")
  print("Distributed Computing & Scalability Test Data:")
  print(runtime_data)
  ggplot(runtime_data, aes(x = Data_Size, y = Runtime, color = Method)) +
    geom_line() +
    geom_point() +
    labs(title = "Runtime Across Different Data Sizes", 
         x = "Data Size", y = "Runtime (seconds)")
  ggsave("runtime_comparison.png")
  write.csv(runtime_data, "E:/CMML/part2/ICA/paper_result/scrublet/runtime_comparison_output.csv", row.names = FALSE)
}

# ========= 5. Trajectory Analysis Using Slingshot =========
# Using the UMAP coordinates and clustering results from the Seurat object,
# infer cell trajectories with Slingshot.
umap_coords <- Embeddings(seurat_obj, "umap")
cluster_labels <- seurat_obj@meta.data$seurat_clusters
doublet_label <- seurat_obj@meta.data$doublet_label  # Assume "doublet" indicates doublet cells

# Color the cells: red for "doublet" and blue for "singlet".
cell_colors <- ifelse(doublet_label == "doublet", "red", "blue")

# Choose the cluster with the most cells as the starting cluster.
init_clust <- names(which.max(table(cluster_labels)))
sling_obj <- slingshot(umap_coords, clusterLabels = cluster_labels, start.clus = init_clust)

# Plot the UMAP with the inferred trajectories.
plot(umap_coords, 
     col = cell_colors,
     pch = 16, asp = 1,
     main = "Trajectory Analysis: Distribution of Doublet Cells (Slingshot)")
curves_list <- slingCurves(sling_obj)
for (i in seq_along(curves_list)) {
  curve_i <- curves_list[[i]]$curve
  if (!is.null(curve_i)) {
    lines(curve_i[,1], curve_i[,2], col = "black", lwd = 2)
  }
}
legend("topright", legend = c("Doublet", "Singlet"), col = c("red", "blue"), pch = 16)

# Optionally, save the trajectory analysis plot.
png(filename = "E:/CMML/part2/ICA/paper_result/scrublet/trajectory_analysis.png", width = 800, height = 600)
plot(umap_coords, 
     col = cell_colors,
     pch = 16, asp = 1,
     main = "Trajectory Analysis: Distribution of Doublet Cells (Slingshot)")
for (i in seq_along(curves_list)) {
  curve_i <- curves_list[[i]]$curve
  if (!is.null(curve_i)) {
    lines(curve_i[,1], curve_i[,2], col = "black", lwd = 2)
  }
}
legend("topright", legend = c("Doublet", "Singlet"), col = c("red", "blue"), pch = 16)
dev.off()

# ========= 6. UMAP Visualization Based on Predicted Scores =========
# Use FeaturePlot to display the predicted doublet score ("predicted_doublet")
# (which should be present in the metadata) on the UMAP plot.
FeaturePlot(seurat_obj, features = "predicted_doublet") +
  ggtitle("UMAP - Predicted Doublet Scores")










library(dplyr)
library(Seurat)
library(ggplot2)

# 1. 读取 SCRUBLET 文件
scrublet_obs <- read.csv("E:/CMML/part2/ICA/paper_result/scrublet/scrublet_obs_results.csv",
                         header = TRUE, row.names = 1, stringsAsFactors = FALSE)
cat("读取的 scrublet_obs (前3行)：\n")
print(head(scrublet_obs, 3))

# 说明：scrublet_obs 包含以下列：
#   nCount_RNA, nFeature_RNA, cell_annotation, scrublet_score, scrublet_doublet_prediction

# ============ 2. 数据类型转换 ============
# 确保 scrublet_score 为数值型
scrublet_obs$scrublet_score <- as.numeric(as.character(scrublet_obs$scrublet_score))
cat("scrublet_score 数值统计：\n")
print(summary(scrublet_obs$scrublet_score))

# ============ 3. 合并到 Seurat 对象 ============
# 假设 seurat_obj 已经基于 count 矩阵创建（细胞名称应与 scrublet_obs 行名一致）
# 检查 Seurat 对象中的细胞名称（例如前6个）
cat("Seurat 对象中的细胞名称 (前6个)：\n")
print(head(Cells(seurat_obj)))

# 检查两个数据集细胞名称的交集
common_cells <- intersect(Cells(seurat_obj), rownames(scrublet_obs))
cat("在 Seurat 对象与 scrublet_obs 中共同存在的细胞数：", length(common_cells), "\n")
if(length(common_cells) == 0){
  stop("合并失败：未找到共同细胞，请检查细胞 ID 是否一致。")
}

# 用 scrublet_obs 文件更新 Seurat 对象的元数据
# 注意：这里假定 seurat_obj 中的细胞顺序与 rownames(scrublet_obs) 匹配
seurat_obj@meta.data <- scrublet_obs[Cells(seurat_obj), ]
cat("合并后的元数据列名：\n")
print(colnames(seurat_obj@meta.data))
print(head(seurat_obj@meta.data, 3))

# ============ 4. 根据 SCRUBLET 得分筛选 Corrected Cells ============
# 计算 60% 分位数作为阈值，降低得分较高（较可能为 doublet）的细胞
threshold_scrub <- quantile(seurat_obj@meta.data$scrublet_score, 0.6, na.rm = TRUE)
cat("SCRUBLET 使用的阈值为：", threshold_scrub, "\n")

# 筛选出预测得分 <= 阈值的细胞（认为为 singlet 的细胞）
filtered_cells <- rownames(seurat_obj@meta.data)[seurat_obj@meta.data$scrublet_score <= threshold_scrub]
cat("筛选出的细胞数：", length(filtered_cells), "\n")

# 检查筛选结果与 Seurat 对象细胞名称的交集
common_cells2 <- intersect(filtered_cells, Cells(seurat_obj))
cat("在 Seurat 对象中匹配到的筛选细胞数：", length(common_cells2), "\n")
if(length(common_cells2) == 0) {
  stop("经过筛选后没有匹配到任何细胞，请检查 scrublet_score 数值与细胞 ID。")
}

# ============ 5. 构建 Corrected Seurat 对象 ============
DefaultAssay(seurat_obj) <- "RNA"  # 确保默认 assay 正确
corrected_seurat <- subset(seurat_obj, cells = common_cells2)
cat("构建的 Corrected Seurat 对象中的细胞数：", ncol(corrected_seurat), "\n")

# ============ 6. 保存 Corrected Seurat 对象的元数据 ============
write.csv(corrected_seurat@meta.data,
          "E:/CMML/part2/ICA/paper_result/scrublet/scrublet_corrected_metadata.csv",
          row.names = TRUE)