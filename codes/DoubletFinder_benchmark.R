# ===== Load Required Packages =====
library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(slingshot)
library(PRROC)
library(data.table)

# ===== 1. Read CSV Files – Parameters, Overall Performance, and Threshold Analysis =====
# (Using DoubletFinder file paths; adjust the file prefixes as needed)
df_obs    <- read.csv("E:/CMML/part2/ICA/paper_result/doubletfinder/doublet_parameters.csv")
df_perf   <- read.csv("E:/CMML/part2/ICA/paper_result/doubletfinder/doublet_overall_performance.csv")
df_thresh <- read.csv("E:/CMML/part2/ICA/paper_result/doubletfinder/doublet_threshold_analysis.csv")

# ===== 2. Differential Expression (DE) Analysis =====
# (Compares doublet versus singlet cells based on true labels stored in "doublet_label")
# Load raw count data (rows: genes, columns: cells) from an RDS file
count_matrix <- readRDS("E:/CMML/part2/ICA/paper_result/counts_matrix.rds")
if (is.null(rownames(count_matrix))) {
  rownames(count_matrix) <- paste0("gene_", 1:nrow(count_matrix))
}
if (is.null(colnames(count_matrix))) {
  colnames(count_matrix) <- paste0("cell_", 1:ncol(count_matrix))
}

# Create a Seurat object using the count matrix and metadata from df_obs
seurat_obj <- CreateSeuratObject(counts = count_matrix, meta.data = df_obs)
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)

# Set the true cell labels (assumed numeric, where 1 indicates a doublet)
labels <- df_obs$doublet_label
labels <- ifelse(labels == 1, "doublet", "singlet")
labels <- as.factor(labels)
seurat_obj@meta.data$doublet_label <- labels
Idents(seurat_obj) <- seurat_obj@meta.data$doublet_label

# Identify differentially expressed genes between doublets and singlets, and save the results
de_results <- FindMarkers(seurat_obj, ident.1 = "doublet", ident.2 = "singlet")
write.csv(de_results, "E:/CMML/part2/ICA/paper_result/doubletfinder/doubletfinder_DE_results.csv", row.names = TRUE)
print("Differential expression analysis completed and saved to doubletfinder_DE_results.csv")

# ===== 3. Clustering and UMAP Visualization =====
# Re-run normalization, variable feature selection, scaling, and PCA if needed
seurat_obj <- NormalizeData(seurat_obj) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(seurat_obj, verbose = FALSE)

# Compute neighbors, clusters, and UMAP (using dimensions 1:10 as an example)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)

# Visualize UMAP with cluster labels
DimPlot(seurat_obj, reduction = "umap", label = TRUE) + ggtitle("Cluster Results")

# Calculate and save the proportion of doublets in each cluster
cluster_doublet_prop <- seurat_obj@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(Doublet_Rate = mean(as.numeric(doublet_label)))
print(cluster_doublet_prop)
write.csv(cluster_doublet_prop, "E:/CMML/part2/ICA/paper_result/doubletfinder/doubletfinder_cluster_doublet_proportions.csv", row.names = FALSE)

# ===== 4. Distributed Computing & Scalability Testing =====
# (If a runtime and memory summary file exists, load and plot the data)
if (file.exists("E:/CMML/part2/ICA/paper_result/runtime_memory_summary.csv")) {
  runtime_data <- fread("E:/CMML/part2/ICA/paper_result/runtime_memory_summary.csv")
  print("Distributed Computing & Scalability Data:")
  print(runtime_data)
  ggplot(runtime_data, aes(x = Data_Size, y = Runtime, color = Method)) +
    geom_line() +
    geom_point() +
    labs(title = "Runtime for Different Data Sizes", x = "Data Size", y = "Runtime (seconds)")
  ggsave("runtime_comparison.png")
  write.csv(runtime_data, "E:/CMML/part2/ICA/paper_result/runtime_comparison_output.csv", row.names = FALSE)
}

# ===== 5. Trajectory Analysis Using Slingshot =====
# Extract UMAP coordinates, cluster labels, and doublet labels from the Seurat object
umap_coords    <- Embeddings(seurat_obj, "umap")
cluster_labels <- seurat_obj@meta.data$seurat_clusters
doublet_label  <- seurat_obj@meta.data$doublet_label

# Set cell colors: "doublet" in red, "singlet" in blue
cell_colors <- ifelse(doublet_label == "doublet", "red", "blue")

# Choose the cluster with the most cells as the initial cluster for trajectory inference
init_clust <- names(which.max(table(cluster_labels)))
sling_obj  <- slingshot(umap_coords, clusterLabels = cluster_labels, start.clus = init_clust)

# Plot the UMAP with cell colors based on true labels
plot(umap_coords, 
     col = cell_colors,
     pch = 16,
     asp = 1,
     main = "Trajectory Analysis: Distribution of Doublet Cells (Slingshot)")
# Overlay the trajectory curves
curves_list <- slingCurves(sling_obj)
for (i in seq_along(curves_list)) {
  curve_i <- curves_list[[i]]$curve
  if (!is.null(curve_i)) {
    lines(curve_i[, 1], curve_i[, 2], col = "black", lwd = 2)
  }
}
legend("topright", legend = c("Doublet", "Singlet"), col = c("red", "blue"), pch = 16)

# Save the trajectory analysis plot as a PNG file
png(filename = "E:/CMML/part2/ICA/paper_result/performance/trajectory_analysis.png", 
    width = 800, height = 600)
plot(umap_coords, 
     col = cell_colors,
     pch = 16,
     asp = 1,
     main = "Trajectory Analysis: Distribution of Doublet Cells (Slingshot)")
for (i in seq_along(curves_list)) {
  curve_i <- curves_list[[i]]$curve
  if (!is.null(curve_i)) {
    lines(curve_i[, 1], curve_i[, 2], col = "black", lwd = 2)
  }
}
legend("topright", legend = c("Doublet", "Singlet"), col = c("red", "blue"), pch = 16)
dev.off()


#
# 设置阈值：以 60% 分位数作为阈值
threshold <- quantile(df_obs$predicted_score, probs = 0.6, na.rm = TRUE)
cat("使用的阈值为：", threshold, "\n")

# 筛选预测得分低于该阈值的细胞
corrected_cells <- rownames(df_obs)[ df_obs$predicted_score <= threshold ]
cat("筛选后细胞数：", length(corrected_cells), "\n")
# 假设 seurat_obj 已利用 doublet_parameters.csv 等数据创建好
# 以下代码将只保留预测为 singlet 的细胞
# 将 corrected_cells 中数字型的字符串转换为数值索引，然后提取对应的细胞名称
corrected_cell_names <- Cells(seurat_obj)[as.numeric(corrected_cells)]
head(corrected_cell_names)
common_cells <- intersect(corrected_cell_names, Cells(seurat_obj))
cat("匹配到的细胞数量：", length(common_cells), "\n")
# 确保默认 assay 为 RNA
DefaultAssay(seurat_obj) <- "RNA"

if (length(common_cells) == 0) {
  stop("没有在 Seurat 对象中找到任何匹配的细胞，请检查细胞名称。")
} else {
  corrected_seurat <- subset(seurat_obj, cells = common_cells)
  cat("构建的 Corrected Seurat 对象中的细胞数量：", ncol(corrected_seurat), "\n")
}
write.csv(corrected_seurat@meta.data, "E:/CMML/part2/ICA/paper_result/doubletfinder/doubletfinder_corrected_metadata.csv", row.names = TRUE)

