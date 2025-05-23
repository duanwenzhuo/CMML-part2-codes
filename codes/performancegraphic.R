# -----------------------------
# Part 1. 加载包、读取数据及数据预处理
# -----------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(PRROC)

# ① 读取真实标签（来自 Solo 的观测结果文件）
solo_obs <- read.csv("E:/CMML/part2/ICA/paper_result/solo/solo_obs_results.csv")
# 假设真实标签存储在 cell_annotation 列中，"doublet" 表示双细胞，其它表示单细胞
solo_obs$true_label <- ifelse(solo_obs$cell_annotation == "doublet", 1, 0)

# ② 读取各方法的预测得分文件

# DoubletFinder：预测得分来自 doublet_parameters.csv
df_df <- read.csv("E:/CMML/part2/ICA/paper_result/doubletfinder/doublet_parameters.csv")
df_df$Method <- "DoubletFinder"

# Scrublet：预测得分所在的列名为 "scrublet_score"
scrublet_df <- read.csv("E:/CMML/part2/ICA/paper_result/scrublet/scrublet_results.csv")
scrublet_df$Method <- "Scrublet"
# 确保将得分转换为数值型
scrublet_df$scrublet_score <- as.numeric(as.character(scrublet_df$scrublet_score))

# Solo：预测得分来自 solo_doublet_scores.csv，注意列名为 "predicted_doublet"
solo_perf <- read.csv("E:/CMML/part2/ICA/paper_result/solo/solo_doublet_scores.csv")
solo_perf$Method <- "Solo"

# ③ 为每个方法赋予真实标签（假设细胞顺序一致，统一使用 solo_obs$true_label）
df_df$true_label <- solo_obs$true_label
scrublet_df$true_label <- solo_obs$true_label
solo_perf$true_label <- solo_obs$true_label

# 同时，确保其它文件中的预测得分为数值型：
df_df$predicted_score <- as.numeric(as.character(df_df$predicted_score))
# 对 Solo，利用 predicted_doublet 列转换为 predicted_score
solo_perf$predicted_score <- as.numeric(as.character(solo_perf$predicted_doublet))

# -----------------------------
# Part 2. 计算 AUROC 与 AUPRC
# -----------------------------
# 对 DoubletFinder
pos_df <- df_df$predicted_score[df_df$true_label == 1]
neg_df <- df_df$predicted_score[df_df$true_label == 0]
pr_df <- pr.curve(scores.class0 = pos_df, scores.class1 = neg_df, curve = TRUE)
roc_df <- roc.curve(scores.class0 = pos_df, scores.class1 = neg_df, curve = TRUE)
AUROC_df <- roc_df$auc
AUPRC_df <- pr_df$auc.integral

# 对 Scrublet（使用 scrublet_score 列）
pos_scrub <- scrublet_df$scrublet_score[scrublet_df$true_label == 1]
neg_scrub <- scrublet_df$scrublet_score[scrublet_df$true_label == 0]
pr_scrub <- pr.curve(scores.class0 = pos_scrub, scores.class1 = neg_scrub, curve = TRUE)
roc_scrub <- roc.curve(scores.class0 = pos_scrub, scores.class1 = neg_scrub, curve = TRUE)
AUROC_scrub <- roc_scrub$auc
AUPRC_scrub <- pr_scrub$auc.integral

# 对 Solo
pos_solo <- solo_perf$predicted_score[solo_perf$true_label == 1]
neg_solo <- solo_perf$predicted_score[solo_perf$true_label == 0]
pr_solo <- pr.curve(scores.class0 = pos_solo, scores.class1 = neg_solo, curve = TRUE)
roc_solo <- roc.curve(scores.class0 = pos_solo, scores.class1 = neg_solo, curve = TRUE)
AUROC_solo <- roc_solo$auc
AUPRC_solo <- pr_solo$auc.integral

# 构建性能数据框（每个方法一行）
performance_df <- data.frame(
  Method = c("DoubletFinder", "Scrublet", "Solo"),
  AUROC = c(AUROC_df, AUROC_scrub, AUROC_solo),
  AUPRC = c(AUPRC_df, AUPRC_scrub, AUPRC_solo)
)
print(performance_df)

# -----------------------------
# Part 3. 绘制 AUROC 与 AUPRC 对比图
# -----------------------------
# 将 performance_df 转换为长格式
perf_long <- pivot_longer(performance_df, cols = c("AUROC", "AUPRC"),
                          names_to = "Metric", values_to = "Value")

p_perf <- ggplot(perf_long, aes(x = Method, y = Value, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~Metric, scales = "free_y") +
  labs(title = "Performance Benchmark: AUROC and AUPRC",
       x = "Method", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_perf)
# 保存图形到指定路径
ggsave("E:/CMML/part2/ICA/paper_result/performance/performance_benchmark.png",
       plot = p_perf, width = 8, height = 6, dpi = 300)

# -----------------------------
# Part 4. 计算不同识别率下的 Precision, Recall 和 TNR
# -----------------------------
# 定义识别率：取预测得分最高的前 10%, 20%, 40% 的细胞
rates <- c(0.1, 0.2, 0.4)

# 对 DoubletFinder 指标计算
metrics_df_df <- data.frame()
for (r in rates) {
  n <- nrow(df_df)
  threshold <- sort(df_df$predicted_score, decreasing = TRUE)[floor(n * r)]
  pred_label <- ifelse(df_df$predicted_score >= threshold, 1, 0)
  TP <- sum(pred_label == 1 & df_df$true_label == 1)
  FP <- sum(pred_label == 1 & df_df$true_label == 0)
  FN <- sum(pred_label == 0 & df_df$true_label == 1)
  TN <- sum(pred_label == 0 & df_df$true_label == 0)
  precision <- ifelse((TP + FP) == 0, NA, TP / (TP + FP))
  recall <- ifelse((TP + FN) == 0, NA, TP / (TP + FN))
  TNR <- ifelse((TN + FP) == 0, NA, TN / (TN + FP))
  metrics_df_df <- rbind(metrics_df_df,
                         data.frame(Method = "DoubletFinder",
                                    IdentificationRate = r * 100,
                                    Precision = precision,
                                    Recall = recall,
                                    TNR = TNR))
}

# 对 Scrublet 指标计算（使用 scrublet_score）
metrics_scrub <- data.frame()
for (r in rates) {
  n <- nrow(scrublet_df)
  threshold <- sort(scrublet_df$scrublet_score, decreasing = TRUE)[floor(n * r)]
  pred_label <- ifelse(scrublet_df$scrublet_score >= threshold, 1, 0)
  TP <- sum(pred_label == 1 & scrublet_df$true_label == 1)
  FP <- sum(pred_label == 1 & scrublet_df$true_label == 0)
  FN <- sum(pred_label == 0 & scrublet_df$true_label == 1)
  TN <- sum(pred_label == 0 & scrublet_df$true_label == 0)
  precision <- ifelse((TP + FP) == 0, NA, TP / (TP + FP))
  recall <- ifelse((TP + FN) == 0, NA, TP / (TP + FN))
  TNR <- ifelse((TN + FP) == 0, NA, TN / (TN + FP))
  metrics_scrub <- rbind(metrics_scrub,
                         data.frame(Method = "Scrublet",
                                    IdentificationRate = r * 100,
                                    Precision = precision,
                                    Recall = recall,
                                    TNR = TNR))
}

# 对 Solo 指标计算
metrics_solo <- data.frame()
for (r in rates) {
  n <- nrow(solo_perf)
  threshold <- sort(solo_perf$predicted_score, decreasing = TRUE)[floor(n * r)]
  pred_label <- ifelse(solo_perf$predicted_score >= threshold, 1, 0)
  TP <- sum(pred_label == 1 & solo_perf$true_label == 1)
  FP <- sum(pred_label == 1 & solo_perf$true_label == 0)
  FN <- sum(pred_label == 0 & solo_perf$true_label == 1)
  TN <- sum(pred_label == 0 & solo_perf$true_label == 0)
  precision <- ifelse((TP + FP) == 0, NA, TP / (TP + FP))
  recall <- ifelse((TP + FN) == 0, NA, TP / (TP + FN))
  TNR <- ifelse((TN + FP) == 0, NA, TN / (TN + FP))
  metrics_solo <- rbind(metrics_solo,
                        data.frame(Method = "Solo",
                                   IdentificationRate = r * 100,
                                   Precision = precision,
                                   Recall = recall,
                                   TNR = TNR))
}

# 合并所有方法的识别率指标
metrics_all <- rbind(metrics_df_df, metrics_scrub, metrics_solo)
print(metrics_all)

# 转换为长格式便于绘图
metrics_long <- pivot_longer(metrics_all, cols = c("Precision", "Recall", "TNR"),
                             names_to = "Metric", values_to = "Value")

p_metrics <- ggplot(metrics_long, aes(x = as.factor(IdentificationRate), y = Value,
                                      group = Method, color = Method)) +
  geom_line(aes(linetype = Method), size = 1.2, na.rm = TRUE) +
  geom_point(size = 3, na.rm = TRUE) +
  facet_wrap(~Metric, scales = "free_y") +
  labs(title = "Metrics under Different Identification Rates",
       x = "Identification Rate (%)", y = "Metric Value") +
  theme_minimal() +
  theme(text = element_text(size = 12))
print(p_metrics)
ggsave("E:/CMML/part2/ICA/paper_result/performance/performance_metrics_identification_rate.png",
       plot = p_metrics, width = 12, height = 4, dpi = 300)
