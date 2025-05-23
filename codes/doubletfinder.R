# 加载必要的包
library(Seurat)
library(DoubletFinder)
library(PRROC)
library(Matrix)
library(ggplot2)
library(data.table)
library(pryr)      # 用于内存检测（模拟分析部分）
library(gridExtra) # 可选：用于绘图组合
# 读取 RDS 文件
data <- readRDS("E:/CMML/part2/ICA/real_datasets/cline-ch.rds")

# 提取计数矩阵和标签
counts <- data[[1]]
labels_raw <- data[[2]]
# 转换标签为0/1形式：doublet为1，singlet为0
labels <- ifelse(labels_raw == "doublet", 1, 0)
# 创建Seurat对象（设置 min.cells 和 min.features 可根据实际情况调整）
seurat_obj <- CreateSeuratObject(counts = counts, min.cells = 3, min.features = 200)

# 数据归一化
seurat_obj <- NormalizeData(seurat_obj)

# 寻找高变基因（vst方法，提取2000个变异基因）
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# 数据标准化
seurat_obj <- ScaleData(seurat_obj)

# PCA计算（可调整PC数量）
seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)

# UMAP降维（使用前10个主成分）
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
# 设置默认assay为RNA
DefaultAssay(seurat_obj) <- "RNA"

# 参数扫描：在1到10个PC上运行（sct设为FALSE）
sweep_res <- paramSweep(seurat_obj, PCs = 1:10, sct = FALSE)
sweep_stats <- summarizeSweep(sweep_res)
optimal_pK <- as.numeric(as.character(find.pK(sweep_stats)$pK[1]))

# 计算预设的双细胞数（nExp），这里根据真实标签计算：双细胞细胞数 = sum(labels)
nExp <- sum(labels)

# 执行DoubletFinder
seurat_obj <- doubletFinder(
  seurat_obj,
  PCs = 1:10,
  pN = 0.25,
  pK = optimal_pK,
  nExp = nExp
)

# 自动获取DoubletFinder计算得分的列名（列名通常以"pANN_"开头）
result_column <- grep("^pANN_", colnames(seurat_obj@meta.data), value = TRUE)
# 提取预测得分
scores <- seurat_obj@meta.data[[result_column]]
# 根据预测得分设定阈值（这里简单假设 >0.5 预测为双细胞）
predicted_labels <- ifelse(scores > 0.5, 1, 0)

# 整理结果数据框
doublet_params <- data.frame(
  doublet_label = predicted_labels,
  predicted_score = scores,
  nExp = nExp
)

# 保存到CSV（请确保路径存在）
write.csv(doublet_params, "E:/CMML/part2/ICA/paper_result/doubletfinder/doublet_parameters.csv", 
          row.names = FALSE)
# 计算前景和背景得分
fg <- scores[labels == 1]
bg <- scores[labels == 0]

# 计算PR和ROC曲线
pr_curve <- pr.curve(scores.class0 = fg, scores.class1 = bg)
roc_curve <- roc.curve(scores.class0 = fg, scores.class1 = bg)

# 保存整体性能指标
performance <- list(
  AUPRC = pr_curve$auc.integral,
  AUROC = roc_curve$auc,
  PR = pr_curve,
  ROC = roc_curve
)

# 对不同识别率（10%, 20%, 40%）进行阈值分析
threshold_rates <- c(0.1, 0.2, 0.4)
threshold_analysis <- sapply(threshold_rates, function(rate) {
  n_cells <- floor(length(labels) * rate)
  threshold_value <- sort(scores, decreasing = TRUE)[n_cells]
  predictions <- scores > threshold_value
  conf_matrix <- table(Predicted = predictions, Actual = labels)
  
  # 计算指标：Precision, Recall, TNR
  tp <- conf_matrix["TRUE", "1"]
  fp <- conf_matrix["TRUE", "0"]
  tn <- conf_matrix["FALSE", "0"]
  fn <- conf_matrix["FALSE", "1"]
  
  precision <- tp / (tp + fp)
  recall <- tp / (tp + fn)
  tnr <- tn / (tn + fp)
  
  c(Precision = precision, Recall = recall, TNR = tnr)
})
# 提取UMAP坐标
umap_coords <- as.data.frame(Embeddings(seurat_obj, "umap"))
umap_coords$cell_id <- rownames(umap_coords)

# 如果Seurat对象元数据中有cell_annotation，则一并导出
if ("cell_annotation" %in% colnames(seurat_obj@meta.data)) {
  umap_coords$doublet_label <- seurat_obj@meta.data$cell_annotation
}

# 将UMAP坐标保存到CSV文件中
write.csv(umap_coords, "E:/CMML/part2/ICA/paper_result/doubletfinder/c", 
          row.names = FALSE)
# 输出性能指标
performance_df <- data.frame(AUROC = performance$AUROC, AUPRC = performance$AUPRC)
write.csv(performance_df, "E:/CMML/part2/ICA/paper_result/doubletfinder/doublet_overall_performance.csv", 
          row.names = FALSE)

# 保存阈值分析结果
threshold_df <- as.data.frame(threshold_analysis)
write.csv(threshold_df, "E:/CMML/part2/ICA/paper_result/doubletfinder/doublet_threshold_analysis.csv", 
          row.names = FALSE)

# 保存整个Seurat对象（后续可用于DE、聚类、轨迹分析等）
saveRDS(seurat_obj, file = "E:/CMML/part2/ICA/paper_result/doubletfinder/seurat_object_with_doublet_detection.rds")

cat("全部结果保存完毕！\n")
# 1. 指定模拟数据规模（细胞数）
num_cells <- 500   # 例如，500个细胞；你可以依次修改为1000, 2000, 4000等
num_genes <- 1000  # 固定基因数

# 2. 生成计数矩阵并设置行列名
count_matrix_sim <- matrix(rpois(num_genes * num_cells, lambda = 10),
                           nrow = num_genes, ncol = num_cells)
rownames(count_matrix_sim) <- paste0("gene_", 1:num_genes)
colnames(count_matrix_sim) <- paste0("cell_", 1:num_cells)

# 转换为稀疏矩阵（确保dimnames不丢失）
sparse_count_matrix <- as(count_matrix_sim, "dgCMatrix")
dimnames(sparse_count_matrix) <- list(rownames(count_matrix_sim), colnames(count_matrix_sim))

# 3. 构造真实标签（10%的细胞为doublet，其余为singlet）
truth <- rep("singlet", num_cells)
doublet_idx <- sample(1:num_cells, size = max(1, floor(0.1 * num_cells)))
truth[doublet_idx] <- "doublet"

# 4. 将模拟数据保存到一个临时RDS文件
temp_file <- tempfile(fileext = ".rds")
saveRDS(list(sparse_count_matrix, truth), temp_file)

# 5. 记录运行时间和内存占用（利用 pryr 包）
time_info <- system.time({
  sim_mem_diff <- pryr::mem_change({
    # 加载模拟数据
    data_sim <- readRDS(temp_file)
    counts_sim <- data_sim[[1]]
    truth_sim <- data_sim[[2]]
    labels_sim <- ifelse(truth_sim == "doublet", 1, 0)
    
    # 构建Seurat对象并跑UMAP
    seurat_obj_sim <- CreateSeuratObject(counts = counts_sim, min.cells = 3, min.features = 200)
    seurat_obj_sim <- NormalizeData(seurat_obj_sim)
    seurat_obj_sim <- FindVariableFeatures(seurat_obj_sim, selection.method = "vst", nfeatures = 2000)
    seurat_obj_sim <- ScaleData(seurat_obj_sim)
    seurat_obj_sim <- RunPCA(seurat_obj_sim, verbose = FALSE)
    seurat_obj_sim <- RunUMAP(seurat_obj_sim, dims = 1:10)
    
    # 双细胞检测
    DefaultAssay(seurat_obj_sim) <- "RNA"
    sweep_res_sim <- paramSweep(seurat_obj_sim, PCs = 1:10, sct = FALSE)
    sweep_stats_sim <- summarizeSweep(sweep_res_sim)
    optimal_pK_sim <- as.numeric(as.character(find.pK(sweep_stats_sim)$pK[1]))
    nExp_sim <- sum(labels_sim)
    seurat_obj_sim <- doubletFinder(seurat_obj_sim, PCs = 1:10, pN = 0.25, pK = optimal_pK_sim, nExp = nExp_sim)
    result_column_sim <- grep("^pANN_", colnames(seurat_obj_sim@meta.data), value = TRUE)
    scores_sim <- seurat_obj_sim@meta.data[[result_column_sim]]
  })
})
mem_diff_mb <- as.numeric(sim_mem_diff) / (1024^2)
runtime_sec <- time_info["elapsed"]

cat("模拟数据 (", num_cells, " cells) 的运行时间：", runtime_sec, "秒，内存增加：", mem_diff_mb, " MB\n")
# -----------------------------------------------
# 追加代码：添加预测双细胞标签并绘制 UMAP 图
# -----------------------------------------------

# -------------------------------
# 追加代码：添加预测标签并绘制UMAP图
# -------------------------------

# 假设您已经完成 DoubletFinder 计算，并在 seurat_obj@meta.data 中获得了得分，
# 我们现在将阈值设为 0.5，生成预测的双细胞标签，并将结果添加到 Seurat 对象中
seurat_obj@meta.data$predicted_doublet_label <- ifelse(seurat_obj@meta.data[[result_column]] > 0.5, "doublet", "singlet")

# 提取UMAP坐标（使用 Embeddings() 函数）
umap_coords <- as.data.frame(Embeddings(seurat_obj, "umap"))
# 添加额外的列：细胞ID和预测的双细胞标签
umap_coords$cell_id <- rownames(umap_coords)
# 这里我们假设之前如果 Seurat 对象中没有 "doublet_label" 列，我们可以用刚刚添加的预测结果替代
umap_coords$doublet_label <- seurat_obj@meta.data$predicted_doublet_label

# 显示列名以确认（保证输出为："umap_1", "umap_2", "cell_id", "doublet_label"）
print(colnames(umap_coords))

# 使用ggplot2绘制UMAP散点图，依据 "umap_1" 和 "umap_2" 绘制，并根据 "doublet_label" 进行上色
library(ggplot2)
p_umap <- ggplot(umap_coords, aes(x = umap_1, y = umap_2, color = doublet_label)) +
  geom_point(alpha = 0.7, size = 1.5) +
  labs(title = "UMAP Plot with DoubletFinder Prediction",
       x = "UMAP1",
       y = "UMAP2",
       color = "Prediction") +
  theme_minimal()

# 显示图形
print(p_umap)

# 保存绘制的图像到指定路径
ggsave("E:/CMML/part2/ICA/paper_result/doubletfinder/umap_plot_with_prediction.png",
       plot = p_umap, width = 8, height = 6, dpi = 300)

# 删除临时文件
unlink(temp_file)
