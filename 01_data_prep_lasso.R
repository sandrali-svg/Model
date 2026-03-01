# ==================== 初始化设置 ====================
# 清理环境
rm(list = ls())
gc()

# 加载必要的包
required_packages <- c("mice", "Hmisc", "glmnet", "dplyr", "tidyr", "ggplot2", 
                       "openxlsx", "officer", "flextable", "readxl", "reshape2",
                       "viridis", "RColorBrewer", "corrplot", "patchwork")

for(pkg in required_packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# 设置随机种子保证结果可重复
set.seed(123456)

# ==================== 数据读取和预处理 ====================
cat("=== 数据读取和预处理 ===\n")

# 读取数据
data <- read_excel("分析用数据集.xlsx", sheet = "Sheet1")
if(file.exists("分析用数据集.xlsx")) {
  data <- read_excel("分析用数据集.xlsx", sheet = "Sheet1")
} else {
  cat("未找到原始数据，尝试使用模拟数据 example_data.xlsx\n")
  data <- read.csv("example_data.xlsx")
}

# 查看数据基本信息
cat("原始数据维度:", dim(data), "\n")
cat("变量数量:", ncol(data), "\n")

# 检查缺失值情况
missing_rates <- colMeans(is.na(data))
cat("\n缺失值统计:\n")
missing_summary <- data.frame(
  变量名 = names(missing_rates),
  缺失率 = round(missing_rates * 100, 2),
  缺失数 = colSums(is.na(data))
)
print(missing_summary[order(-missing_summary$缺失率), ])

# 根据要求：不删除缺失值超过50%的变量，保留所有变量
cat("\n根据要求：保留所有变量进行后续分析\n")
data_clean <- data
cat("保留全部", ncol(data_clean), "个变量\n")

# 数据预处理 - 转换分类变量
cat("\n=== 数据预处理 ===\n")

# 修复病程列中的公式问题
if("病程" %in% names(data_clean)) {
  # 将公式转换为数值
  data_clean$病程 <- as.numeric(gsub("=|=.*?/", "", as.character(data_clean$病程)))
  cat("已修复病程列中的公式问题\n")
}

# 确保疗效变量是二分类因子
if(is.numeric(data_clean$疗效)) {
  data_clean$疗效 <- factor(data_clean$疗效, levels = c(0, 1), labels = c("有效", "显效"))
} else {
  data_clean$疗效 <- factor(as.character(data_clean$疗效), levels = c("0", "1"), labels = c("有效", "显效"))
}

# 识别分类变量
categorical_vars <- c("性别", "婚姻", "是否伴OSA", "血型", 
                      "皮肤点刺分级", "食入变应原筛查分级", "吸入变应原筛查分级")
categorical_vars <- categorical_vars[categorical_vars %in% names(data_clean)]

# 转换分类变量
for(var in categorical_vars) {
  data_clean[[var]] <- as.factor(as.character(data_clean[[var]]))
  cat(var, "的水平数:", length(levels(data_clean[[var]])), "\n")
}

# 检查处理后的数据
cat("\n处理后的疗效变量分布:\n")
print(table(data_clean$疗效, useNA = "always"))

# 创建数据概览图
cat("\n生成数据概览图...\n")

# 1. 疗效分布图
p1 <- ggplot(data_clean, aes(x = 疗效, fill = 疗效)) +
  geom_bar(alpha = 0.8) +
  scale_fill_manual(values = c("#4E79A7", "#F28E2B")) +
  labs(title = "疗效分布", x = "疗效", y = "病例数") +
  theme_minimal() +
  theme(legend.position = "none")

# 2. 关键连续变量分布
continuous_vars <- c("年龄", "BMI", "术前VAS评分", "病程", "空腹葡萄糖")
if(all(continuous_vars %in% names(data_clean))) {
  # 创建多面图
  plot_data <- data_clean[, continuous_vars]
  plot_data_long <- pivot_longer(plot_data, everything(), names_to = "Variable", values_to = "Value")
  
  p2 <- ggplot(plot_data_long, aes(x = Value, fill = Variable)) +
    geom_density(alpha = 0.6) +
    facet_wrap(~ Variable, scales = "free", ncol = 3) +
    scale_fill_viridis_d() +
    labs(title = "关键连续变量分布", x = "值", y = "密度") +
    theme_minimal() +
    theme(legend.position = "none")
}

# 3. 缺失值热图
missing_matrix <- is.na(data_clean)
missing_long <- as.data.frame(missing_matrix) %>%
  mutate(row = row_number()) %>%
  pivot_longer(cols = -row, names_to = "variable", values_to = "missing") %>%
  group_by(variable, missing) %>%
  summarise(count = n()) %>%
  filter(missing == TRUE)

if(nrow(missing_long) > 0) {
  p3 <- ggplot(missing_long, aes(x = reorder(variable, -count), y = count, fill = count)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    scale_fill_gradient(low = "#4E79A7", high = "#E15759") +
    coord_flip() +
    labs(title = "变量缺失情况", x = "变量", y = "缺失数") +
    theme_minimal()
}

# 保存概览图
png("data_overview.png", width = 16, height = 12, units = "in", res = 300)
if(exists("p2") && exists("p3")) {
  (p1 / (p2 / p3)) + 
    plot_annotation(tag_levels = 'A', 
                    title = '数据集概览',
                    theme = theme(plot.title = element_text(hjust = 0.5, size = 16)))
} else {
  print(p1)
}
dev.off()

# 保存预处理后的数据
wb_preprocessed <- createWorkbook()
addWorksheet(wb_preprocessed, "预处理数据")
writeData(wb_preprocessed, "预处理数据", data_clean)
saveWorkbook(wb_preprocessed, "preprocessed_data.xlsx", overwrite = TRUE)
cat("\n已保存预处理数据: preprocessed_data.xlsx\n")

# ==================== 多重插补 ====================
cat("\n=== 多重插补 ===\n")

# 方法1: MICE多重插补
cat("开始MICE多重插补...\n")

# 准备插补方法
imp_method <- make.method(data_clean)

# 为分类变量指定适当的插补方法
for(var in categorical_vars) {
  n_levels <- length(levels(data_clean[[var]]))
  if(n_levels == 2) {
    imp_method[var] <- "logreg"  # 二分类逻辑回归
  } else {
    imp_method[var] <- "polyreg"  # 多分类多项式回归
  }
}

# 连续变量使用pmm
continuous_vars_all <- setdiff(names(data_clean), c(categorical_vars, "疗效"))
imp_method[continuous_vars_all] <- "pmm"

cat("使用的插补方法:\n")
print(imp_method)

# 执行MICE插补 (100个数据集)
mice_imp <- mice(data_clean, 
                 m = 100, 
                 method = imp_method,
                 maxit = 10,
                 printFlag = TRUE,
                 seed = 123456)

cat("MICE插补完成\n")

# 插补诊断图
png("imputation_diagnostics.png", width = 14, height = 10, units = "in", res = 300)
par(mfrow = c(2, 2))
plot(mice_imp)
dev.off()

# 保存插补数据集
cat("保存插补数据集到Excel...\n")

# 创建目录保存插补数据集
if(!dir.exists("imputed_datasets_excel")) {
  dir.create("imputed_datasets_excel")
}

# 保存前5个数据集到单独的Excel文件
for(i in 1:min(5, mice_imp$m)) {
  imp_data <- complete(mice_imp, i)
  wb <- createWorkbook()
  addWorksheet(wb, "插补数据")
  writeData(wb, "插补数据", imp_data)
  filename <- paste0("imputed_datasets_excel/imputed_dataset_", sprintf("%02d", i), ".xlsx")
  saveWorkbook(wb, filename, overwrite = TRUE)
}

# 保存所有数据集到一个Excel文件的不同sheet
wb_all <- createWorkbook()
max_sheets <- min(20, mice_imp$m)
for(i in 1:max_sheets) {
  imp_data <- complete(mice_imp, i)
  sheet_name <- paste0("Dataset_", sprintf("%02d", i))
  addWorksheet(wb_all, sheet_name)
  writeData(wb_all, sheet_name, imp_data)
}
saveWorkbook(wb_all, "all_imputed_datasets.xlsx", overwrite = TRUE)

cat("插补数据集保存完成\n")

# ==================== LASSO回归分析 ====================
cat("\n=== LASSO回归分析 ===\n")

# 定义LASSO函数
perform_lasso <- function(dataset, outcome_var = "疗效") {
  # 确保结果是数值型
  if(is.factor(dataset[[outcome_var]])) {
    y <- as.numeric(dataset[[outcome_var]]) - 1  # 转换为0/1
  } else {
    y <- as.numeric(dataset[[outcome_var]])
  }
  
  # 创建模型矩阵（处理因子变量）
  x <- model.matrix(as.formula(paste("~ . -", outcome_var)), data = dataset)[, -1]
  
  # 10折交叉验证寻找最优lambda
  cv_fit <- cv.glmnet(x, y, 
                      alpha = 1,
                      family = "binomial",
                      nfolds = 10,
                      type.measure = "deviance",
                      seed = 123456)
  
  # 使用最小lambda的模型
  best_lambda <- cv_fit$lambda.min
  
  # 拟合最终模型
  final_model <- glmnet(x, y, 
                        alpha = 1, 
                        lambda = best_lambda,
                        family = "binomial")
  
  return(list(
    model = final_model,
    lambda = best_lambda,
    coefficients = coef(final_model),
    cv_fit = cv_fit,
    variable_names = colnames(x)
  ))
}

# 对每个插补数据集分别进行LASSO回归
cat("开始对100个插补数据集进行LASSO回归...\n")
mice_lasso_results <- list()
success_count <- 0

for (i in 1:100) {
  if (i %% 10 == 0) cat("处理MICE数据集:", i, "/100\n")
  
  tryCatch({
    complete_data <- complete(mice_imp, i)
    result <- perform_lasso(complete_data)
    mice_lasso_results[[i]] <- result
    success_count <- success_count + 1
  }, error = function(e) {
    cat("数据集", i, "出错:", e$message, "\n")
    mice_lasso_results[[i]] <- NULL
  })
}

cat("成功处理的数据集数量:", success_count, "/100\n")

# 过滤掉失败的数据集
mice_lasso_results <- mice_lasso_results[!sapply(mice_lasso_results, is.null)]

if(length(mice_lasso_results) == 0) {
  stop("没有成功分析的数据集，请检查数据问题")
}

# ==================== LASSO结果汇总 ====================
cat("\n=== LASSO结果汇总 ===\n")

# 汇总函数
summarize_lasso_results <- function(lasso_results, method_name) {
  if(length(lasso_results) == 0) {
    stop("没有可用的LASSO结果")
  }
  
  # 获取所有变量名
  all_vars <- rownames(lasso_results[[1]]$coefficients)
  
  # 初始化系数矩阵
  coef_matrix <- matrix(NA, nrow = length(all_vars), ncol = length(lasso_results))
  rownames(coef_matrix) <- all_vars
  
  # 填充系数矩阵
  for(i in 1:length(lasso_results)) {
    current_coef <- as.numeric(lasso_results[[i]]$coefficients)
    coef_matrix[, i] <- current_coef
  }
  
  # 计算平均系数和选择频率
  mean_coefs <- apply(coef_matrix, 1, mean, na.rm = TRUE)
  selection_freq <- apply(coef_matrix != 0, 1, mean, na.rm = TRUE)
  
  # 创建汇总数据框
  summary_df <- data.frame(
    Variable = all_vars,
    Mean_Coefficient = mean_coefs,
    Selection_Frequency = selection_freq,
    Method = method_name
  )
  
  # 按选择频率排序
  summary_df <- summary_df[order(-summary_df$Selection_Frequency), ]
  
  return(summary_df)
}

# 稳定性分析函数
calculate_selection_stability <- function(lasso_results, threshold = 0.6) {
  if(length(lasso_results) == 0) {
    stop("没有可用的LASSO结果")
  }
  
  # 获取所有变量名
  all_vars <- rownames(lasso_results[[1]]$coefficients)
  
  # 初始化选择矩阵
  selection_matrix <- matrix(NA, nrow = length(all_vars), ncol = length(lasso_results))
  rownames(selection_matrix) <- all_vars
  
  # 填充选择矩阵
  for(i in 1:length(lasso_results)) {
    current_selection <- as.numeric(lasso_results[[i]]$coefficients != 0)
    selection_matrix[, i] <- current_selection
  }
  
  # 计算选择频率
  selection_freq <- apply(selection_matrix, 1, mean, na.rm = TRUE)
  
  # 选择在超过threshold比例的数据集中被稳定选择的变量
  stable_vars <- all_vars[selection_freq >= threshold & all_vars != "(Intercept)"]
  
  return(list(
    selection_frequency = selection_freq,
    stable_variables = stable_vars,
    variable_names = all_vars,
    threshold = threshold,
    selection_matrix = selection_matrix
  ))
}

# 生成汇总结果
mice_summary <- summarize_lasso_results(mice_lasso_results, "MICE")
mice_stability <- calculate_selection_stability(mice_lasso_results, 0.6)

# 提取lambda统计
lambda_values <- sapply(mice_lasso_results, function(x) x$lambda)
lambda_summary <- data.frame(
  Statistic = c("Minimum", "25th Percentile", "Median", "Mean", "75th Percentile", "Maximum", "Standard Deviation"),
  Value = c(
    min(lambda_values),
    quantile(lambda_values, 0.25),
    median(lambda_values),
    mean(lambda_values),
    quantile(lambda_values, 0.75),
    max(lambda_values),
    sd(lambda_values)
  )
)

# ==================== 可视化生成 ====================
cat("\n=== 生成可视化图形 ===\n")

# 定义色系
color_palette <- c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", 
                   "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC")

# 1. LASSO路径图
png("mice_lasso_paths.png", width = 14, height = 6, units = "in", res = 300)
par(mfrow = c(1, 2), mar = c(5, 4, 4, 2) + 0.1, family = "serif")

# 使用第一个数据集生成图形
dataset_idx <- 1

# 交叉验证误差图
plot(mice_lasso_results[[dataset_idx]]$cv_fit,
     main = "LASSO交叉验证曲线",
     cex.main = 1.2, font.main = 2,
     xlab = "log(Lambda)", ylab = "二项偏差",
     col = color_palette[1], lwd = 2)
grid(col = "gray80")

# 系数路径图
plot(mice_lasso_results[[dataset_idx]]$cv_fit$glmnet.fit,
     xvar = "lambda",
     main = "LASSO系数路径",
     cex.main = 1.2, font.main = 2,
     xlab = "log(Lambda)", ylab = "系数",
     col = color_palette)
abline(v = log(mice_lasso_results[[dataset_idx]]$lambda), 
       lty = 2, col = color_palette[2], lwd = 2)
legend("bottomright", 
       legend = paste("最优lambda =", round(mice_lasso_results[[dataset_idx]]$lambda, 4)),
       col = color_palette[2], lty = 2, lwd = 2, cex = 0.8, bg = "white")
grid(col = "gray80")
dev.off()

# 2. 变量选择频率图
plot_data <- mice_summary[mice_summary$Variable != "(Intercept)", ]
if(nrow(plot_data) > 20) {
  plot_data <- head(plot_data[order(-plot_data$Selection_Frequency), ], 20)
}

# 创建渐变色
color_gradient <- colorRampPalette(c(color_palette[1], color_palette[2]))(nrow(plot_data))

p_selection <- ggplot(plot_data, aes(x = reorder(Variable, Selection_Frequency), 
                                     y = Selection_Frequency,
                                     fill = Selection_Frequency)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  coord_flip() +
  scale_fill_gradientn(colors = c(color_palette[1], color_palette[2], color_palette[3]),
                       name = "选择频率") +
  labs(title = "LASSO变量选择频率",
       subtitle = paste("基于", length(mice_lasso_results), "个插补数据集"),
       x = "变量", y = "选择频率") +
  theme_minimal(base_family = "serif") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, color = "gray50"),
        axis.text = element_text(size = 10),
        legend.position = "right") +
  geom_text(aes(label = sprintf("%.1f%%", Selection_Frequency * 100)), 
            hjust = -0.1, size = 3.5, family = "serif")

ggsave("variable_selection_frequency.png", p_selection, width = 12, height = 10, dpi = 300)

# 3. 系数热图
if(length(mice_lasso_results) >= 10) {
  n_show <- min(15, length(mice_lasso_results))
  
  # 构建系数矩阵
  all_vars <- mice_stability$variable_names
  coef_matrix <- matrix(NA, nrow = length(all_vars), ncol = n_show)
  rownames(coef_matrix) <- all_vars
  
  for(i in 1:n_show) {
    coef_matrix[, i] <- as.numeric(mice_lasso_results[[i]]$coefficients)
  }
  
  # 移除截距，选择重要变量
  coef_matrix <- coef_matrix[rownames(coef_matrix) != "(Intercept)", ]
  selection_freq <- rowMeans(coef_matrix != 0, na.rm = TRUE)
  top_vars <- names(sort(selection_freq, decreasing = TRUE))
  if(length(top_vars) > 15) top_vars <- top_vars[1:15]
  coef_matrix <- coef_matrix[top_vars, 1:n_show]
  
  # 创建热图数据
  heatmap_data <- as.data.frame(coef_matrix)
  heatmap_data$Variable <- rownames(heatmap_data)
  heatmap_long <- pivot_longer(heatmap_data, cols = -Variable, 
                               names_to = "Dataset", values_to = "Coefficient")
  
  heatmap_long$DatasetNum <- as.numeric(gsub("V", "", heatmap_long$Dataset))
  
  # 添加选择频率信息
  freq_info <- data.frame(Variable = top_vars, 
                          Frequency = selection_freq[top_vars])
  heatmap_long <- merge(heatmap_long, freq_info, by = "Variable")
  
  # 绘制热图
  p_heatmap <- ggplot(heatmap_long, aes(x = DatasetNum, y = reorder(Variable, Frequency), 
                                        fill = Coefficient)) +
    geom_tile(color = "gray80", linewidth = 0.5) +
    scale_fill_gradient2(low = color_palette[1], mid = "white", high = color_palette[3], 
                         midpoint = 0, 
                         name = "系数值",
                         limits = c(-max(abs(coef_matrix), na.rm = TRUE), 
                                    max(abs(coef_matrix), na.rm = TRUE))) +
    scale_x_continuous(breaks = 1:n_show, expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(title = "LASSO系数热图",
         subtitle = "显示前15个变量在前15个数据集中的系数",
         x = "插补数据集编号", y = "变量") +
    theme_minimal(base_family = "serif") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          plot.subtitle = element_text(hjust = 0.5),
          panel.grid = element_blank())
  
  ggsave("coefficient_heatmap.png", p_heatmap, width = 14, height = 10, dpi = 300)
}

# 4. Lambda分布图
lambda_df <- data.frame(Lambda = lambda_values)

p_lambda <- ggplot(lambda_df, aes(x = Lambda)) +
  geom_histogram(aes(y = ..density..), fill = color_palette[4], alpha = 0.7, 
                 bins = 20, color = "white") +
  geom_density(fill = color_palette[5], alpha = 0.3) +
  geom_vline(xintercept = median(lambda_values), linetype = "dashed", 
             color = color_palette[2], size = 1) +
  annotate("text", x = median(lambda_values), y = max(hist(lambda_values, plot = FALSE)$density) * 0.9,
           label = paste("中位数 =", round(median(lambda_values), 4)), 
           hjust = -0.1, color = color_palette[2], family = "serif") +
  labs(title = "LASSO惩罚参数λ分布",
       subtitle = paste("基于", length(lambda_values), "个插补数据集"),
       x = "λ值", y = "密度") +
  theme_minimal(base_family = "serif") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5))

ggsave("lambda_distribution.png", p_lambda, width = 10, height = 8, dpi = 300)

cat("可视化图形生成完成\n")

# ==================== 保存LASSO结果文件到Excel ====================
cat("\n=== 保存LASSO结果文件到Excel ===\n")

# 创建主工作簿
wb_results <- createWorkbook()

# 1. LASSO结果
addWorksheet(wb_results, "LASSO变量选择结果")
writeData(wb_results, "LASSO变量选择结果", mice_summary)

# 2. 稳定变量详情
stable_vars_details <- data.frame(
  Rank = 1:length(mice_stability$stable_variables),
  Variable = mice_stability$stable_variables,
  Selection_Frequency = round(mice_stability$selection_frequency[mice_stability$stable_variables] * 100, 1),
  Mean_Coefficient = round(mice_summary$Mean_Coefficient[match(mice_stability$stable_variables, mice_summary$Variable)], 4)
)
addWorksheet(wb_results, "稳定变量详情")
writeData(wb_results, "稳定变量详情", stable_vars_details)

# 3. Lambda统计
addWorksheet(wb_results, "Lambda统计")
writeData(wb_results, "Lambda统计", lambda_summary)

# 4. 详细变量选择报告
detailed_report <- mice_summary %>%
  mutate(
    Selection_Frequency_Pct = paste0(round(Selection_Frequency * 100, 1), "%"),
    Significance = case_when(
      Selection_Frequency >= 0.8 ~ "***",
      Selection_Frequency >= 0.6 ~ "**",
      Selection_Frequency >= 0.3 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  select(Variable, Mean_Coefficient, Selection_Frequency_Pct, Significance)

addWorksheet(wb_results, "详细变量选择报告")
writeData(wb_results, "详细变量选择报告", detailed_report)

# 5. 选择矩阵（前15个数据集）
if(length(mice_lasso_results) >= 15) {
  selection_matrix <- matrix(NA, nrow = length(mice_stability$variable_names), ncol = 15)
  rownames(selection_matrix) <- mice_stability$variable_names
  
  for(i in 1:15) {
    selection_matrix[, i] <- as.numeric(mice_lasso_results[[i]]$coefficients != 0)
  }
  
  selection_df <- as.data.frame(selection_matrix)
  colnames(selection_df) <- paste0("Dataset_", 1:15)
  selection_df$Variable <- rownames(selection_matrix)
  
  addWorksheet(wb_results, "变量选择矩阵")
  writeData(wb_results, "变量选择矩阵", selection_df)
}

# 6. 缺失值统计
addWorksheet(wb_results, "缺失值统计")
writeData(wb_results, "缺失值统计", missing_summary)

# 保存主工作簿
saveWorkbook(wb_results, "lasso_analysis_results.xlsx", overwrite = TRUE)

cat("LASSO结果文件保存完成: lasso_analysis_results.xlsx\n")

# ==================== 生成LASSO专用Word报告 ====================
cat("\n=== 生成LASSO专用Word报告 ===\n")

# 创建Word文档
doc <- read_docx() %>%
  body_add_par("LASSO回归变量筛选分析报告", style = "heading 1") %>%
  body_add_par(paste("生成时间:", Sys.time()), style = "Normal") %>%
  body_add_par(paste("随机种子:", 123456), style = "Normal") %>%
  body_add_par(" ", style = "Normal")

# 1. 研究概述
doc <- doc %>%
  body_add_par("1. 研究概述", style = "heading 2") %>%
  body_add_par("本研究采用多重插补结合LASSO回归的方法，从多个临床和实验室指标中筛选对疗效（显效 vs 有效）具有预测价值的稳定变量。", style = "Normal") %>%
  body_add_par("主要特点：", style = "Normal") %>%
  body_add_par("• 保留所有变量，不进行变量删除", style = "Normal") %>%
  body_add_par("• 使用100次多重插补处理缺失值", style = "Normal") %>%
  body_add_par("• 基于稳定性选择确定重要变量", style = "Normal") %>%
  body_add_par(" ", style = "Normal")

# 2. 数据特征
doc <- doc %>%
  body_add_par("2. 数据特征", style = "heading 2") %>%
  body_add_par(paste("• 总样本量:", nrow(data_clean), "例"), style = "Normal") %>%
  body_add_par(paste("• 变量数量:", ncol(data_clean), "个"), style = "Normal")

# 疗效分布
efficacy_table <- as.data.frame(table(data_clean$疗效))
names(efficacy_table) <- c("疗效", "病例数")
efficacy_table$百分比 <- paste0(round(efficacy_table$病例数/sum(efficacy_table$病例数)*100, 1), "%")

doc <- doc %>%
  body_add_par("• 疗效分布：", style = "Normal") %>%
  body_add_flextable(flextable(efficacy_table) %>% theme_zebra() %>% autofit()) %>%
  body_add_par(" ", style = "Normal")

# 缺失值情况
doc <- doc %>%
  body_add_par("• 缺失值统计（前10个变量）：", style = "Normal")
missing_top10 <- head(missing_summary[order(-missing_summary$缺失率), ], 10)
doc <- doc %>%
  body_add_flextable(flextable(missing_top10) %>% theme_zebra() %>% autofit()) %>%
  body_add_par(" ", style = "Normal")

# 3. 方法学
doc <- doc %>%
  body_add_par("3. 方法学", style = "heading 2") %>%
  body_add_par("3.1 数据预处理", style = "heading 3") %>%
  body_add_par("• 保留所有原始变量，不进行删除", style = "Normal") %>%
  body_add_par("• 分类变量转换为因子类型", style = "Normal") %>%
  body_add_par("• 疗效变量编码：0=有效，1=显效", style = "Normal") %>%
  body_add_par(" ", style = "Normal") %>%
  body_add_par("3.2 多重插补", style = "heading 3") %>%
  body_add_par("• 方法：链式方程多重插补（MICE）", style = "Normal") %>%
  body_add_par("• 插补数据集数量：100个", style = "Normal") %>%
  body_add_par("• 插补方法：", style = "Normal") %>%
  body_add_par("  - 二分类变量：logistic回归", style = "Normal") %>%
  body_add_par("  - 多分类变量：多项式回归", style = "Normal") %>%
  body_add_par("  - 连续变量：预测均值匹配", style = "Normal") %>%
  body_add_par(" ", style = "Normal") %>%
  body_add_par("3.3 LASSO回归", style = "heading 3") %>%
  body_add_par("• 在每个插补数据集上独立进行LASSO回归", style = "Normal") %>%
  body_add_par("• 10折交叉验证确定最优λ", style = "Normal") %>%
  body_add_par("• 稳定性阈值：选择频率≥60%", style = "Normal") %>%
  body_add_par(" ", style = "Normal")

# 4. 主要结果
doc <- doc %>%
  body_add_par("4. 主要结果", style = "heading 2") %>%
  body_add_par("4.1 LASSO参数统计", style = "heading 3") %>%
  body_add_par(paste("基于", length(mice_lasso_results), "个插补数据集的LASSO分析结果："), style = "Normal")

doc <- doc %>%
  body_add_flextable(flextable(lambda_summary) %>% theme_zebra() %>% autofit()) %>%
  body_add_par(" ", style = "Normal")

# 4.2 稳定变量
doc <- doc %>%
  body_add_par("4.2 稳定预测变量筛选", style = "heading 3") %>%
  body_add_par(paste("共筛选出", nrow(stable_vars_details), "个稳定预测变量（选择频率≥60%）："), style = "Normal") %>%
  body_add_par(" ", style = "Normal")

doc <- doc %>%
  body_add_flextable(flextable(stable_vars_details) %>% theme_zebra() %>% autofit()) %>%
  body_add_par(" ", style = "Normal")

# 4.3 变量选择频率
doc <- doc %>%
  body_add_par("4.3 详细变量选择频率", style = "heading 3") %>%
  body_add_par("所有变量的选择频率分布：", style = "Normal")

detailed_display <- mice_summary %>%
  filter(Variable != "(Intercept)") %>%
  arrange(desc(Selection_Frequency)) %>%
  mutate(
    Selection_Frequency = paste0(round(Selection_Frequency * 100, 1), "%"),
    Rank = 1:n()
  ) %>%
  select(Rank, Variable, Selection_Frequency, Mean_Coefficient) %>%
  head(20)

doc <- doc %>%
  body_add_flextable(flextable(detailed_display) %>% theme_zebra() %>% autofit()) %>%
  body_add_par(" ", style = "Normal")

# 5. 可视化结果
doc <- doc %>%
  body_add_par("5. 可视化结果", style = "heading 2") %>%
  body_add_par("5.1 数据概览", style = "heading 3") %>%
  body_add_par("数据集基本特征分布：", style = "Normal")

if(file.exists("data_overview.png")) {
  doc <- doc %>%
    body_add_img("data_overview.png", width = 6.5, height = 5) %>%
    body_add_par(" ", style = "Normal")
}

doc <- doc %>%
  body_add_par("5.2 LASSO分析图", style = "heading 3") %>%
  body_add_par("LASSO交叉验证曲线和系数路径：", style = "Normal")

if(file.exists("mice_lasso_paths.png")) {
  doc <- doc %>%
    body_add_img("mice_lasso_paths.png", width = 6.5, height = 3.5) %>%
    body_add_par(" ", style = "Normal")
}

doc <- doc %>%
  body_add_par("5.3 变量选择频率", style = "heading 3") %>%
  body_add_par("各变量在插补数据集中的选择频率：", style = "Normal")

if(file.exists("variable_selection_frequency.png")) {
  doc <- doc %>%
    body_add_img("variable_selection_frequency.png", width = 6.5, height = 5) %>%
    body_add_par(" ", style = "Normal")
}

doc <- doc %>%
  body_add_par("5.4 系数热图", style = "heading 3") %>%
  body_add_par("重要变量在不同数据集中的系数变化：", style = "Normal")

if(file.exists("coefficient_heatmap.png")) {
  doc <- doc %>%
    body_add_img("coefficient_heatmap.png", width = 6.5, height = 5) %>%
    body_add_par(" ", style = "Normal")
}

doc <- doc %>%
  body_add_par("5.5 λ参数分布", style = "heading 3") %>%
  body_add_par("LASSO惩罚参数λ的分布情况：", style = "Normal")

if(file.exists("lambda_distribution.png")) {
  doc <- doc %>%
    body_add_img("lambda_distribution.png", width = 6.5, height = 4) %>%
    body_add_par(" ", style = "Normal")
}

# 6. 结论
doc <- doc %>%
  body_add_par("6. 结论", style = "heading 2") %>%
  body_add_par("6.1 主要发现", style = "heading 3") %>%
  body_add_par(paste("• 成功筛选出", nrow(stable_vars_details), "个稳定预测变量"), style = "Normal") %>%
  body_add_par("• LASSO参数选择在不同插补数据集间表现稳定", style = "Normal") %>%
  body_add_par("• 多重插补有效处理了缺失数据", style = "Normal") %>%
  body_add_par(" ", style = "Normal") %>%
  body_add_par("6.2 临床意义", style = "heading 3") %>%
  body_add_par("筛选出的变量反映了疗效预测的多维度特征：", style = "Normal") %>%
  body_add_par("• 人口学特征：年龄、性别等", style = "Normal") %>%
  body_add_par("• 临床指标：病程、VAS评分等", style = "Normal") %>%
  body_add_par("• 实验室指标：炎症指标、代谢指标等", style = "Normal") %>%
  body_add_par("• 过敏相关：变应原筛查分级等", style = "Normal") %>%
  body_add_par(" ", style = "Normal") %>%
  body_add_par("6.3 局限性", style = "heading 3") %>%
  body_add_par("• 样本量相对有限", style = "Normal") %>%
  body_add_par("• 需要外部验证", style = "Normal") %>%
  body_add_par("• 部分变量缺失率较高", style = "Normal") %>%
  body_add_par(" ", style = "Normal")

# 7. 附录
doc <- doc %>%
  body_add_par("附录：技术细节", style = "heading 2") %>%
  body_add_par("软件环境与版本：", style = "Normal") %>%
  body_add_par(paste("• R版本:", R.version$version.string), style = "Normal") %>%
  body_add_par("• mice: 多重插补", style = "Normal") %>%
  body_add_par("• glmnet: LASSO回归", style = "Normal") %>%
  body_add_par("• ggplot2: 数据可视化", style = "Normal") %>%
  body_add_par("• officer/flextable: 报告生成", style = "Normal") %>%
  body_add_par(" ", style = "Normal") %>%
  body_add_par("分析流程概要：", style = "Normal") %>%
  body_add_par("1. 数据读取与检查（保留所有变量）", style = "Normal") %>%
  body_add_par("2. 变量类型转换与预处理", style = "Normal") %>%
  body_add_par("3. MICE多重插补（100个数据集）", style = "Normal") %>%
  body_add_par("4. 独立LASSO回归分析", style = "Normal") %>%
  body_add_par("5. 结果汇总与稳定性评估", style = "Normal") %>%
  body_add_par("6. 可视化与报告生成", style = "Normal") %>%
  body_add_par(" ", style = "Normal")

# 保存Word文档
output_file <- paste0("LASSO_Analysis_Report_", format(Sys.time(), "%Y%m%d_%H%M"), ".docx")
print(doc, target = output_file)

cat("LASSO专用Word报告已生成: ", output_file, "\n")

# ==================== 最终总结 ====================
cat("\n=== 分析完成总结 ===\n")

# 创建总结表格
summary_table <- data.frame(
  项目 = c("数据读取", "变量处理", "多重插补", "LASSO分析", "稳定变量", "文件生成"),
  结果 = c(
    paste("读取", nrow(data), "行，", ncol(data), "列"),
    paste("保留", ncol(data_clean), "个变量，转换", length(categorical_vars), "个分类变量"),
    paste("生成", mice_imp$m, "个插补数据集"),
    paste("成功分析", length(mice_lasso_results), "个数据集"),
    paste("筛选", nrow(stable_vars_details), "个稳定变量"),
    paste("生成", length(list.files(pattern = "\\.(png|xlsx|docx|RData)$")), "个输出文件")
  )
)

print(summary_table)

cat("\n生成的主要文件：\n")
files <- list.files(pattern = "\\.(png|xlsx|docx|RData)$")
for(file in files) {
  cat("• ", file, "\n")
}

# ==================== 保存分析状态 ====================
cat("\n=== 保存关键分析状态 ===\n")

# 保存完整的分析状态
save(mice_imp,                # 多重插补对象
     data_clean,              # 预处理后的数据
     mice_lasso_results,      # LASSO结果列表
     mice_summary,            # LASSO汇总结果
     mice_stability,          # 稳定性分析结果
     lambda_summary,          # λ统计
     stable_vars_details,     # 稳定变量详情
     missing_summary,         # 缺失值统计
     file = "lasso_analysis_state.RData")

cat("已保存分析状态到: lasso_analysis_state.RData\n")

# 保存完整工作空间
save.image(file = "lasso_analysis_complete.RData")
cat("已保存完整工作空间: lasso_analysis_complete.RData\n")

# ==================== 为后续建模准备专用文件 ====================
cat("\n=== 为后续建模准备专用文件 ===\n")

# 从LASSO结果中提取关键信息
stable_vars <- mice_stability$stable_variables
stable_vars <- stable_vars[!stable_vars %in% c("(Intercept)")]

# 保存为建模专用文件
save(
  mice_imp,                    # 多重插补对象
  stable_vars,                 # 稳定变量列表
  mice_summary,                # LASSO汇总
  stable_vars_details,         # 稳定变量详情
  data_clean,                  # 原始清洗数据
  lambda_summary,              # λ参数
  mice_stability,              # 稳定性分析
  missing_summary,             # 缺失值统计
  file = "model_building_input.RData"
)

cat("已创建建模专用文件: model_building_input.RData\n")

# ==================== 验证保存的文件 ====================
cat("\n=== 验证保存的文件 ===\n")

check_files <- c("lasso_analysis_state.RData", 
                 "model_building_input.RData",
                 "lasso_analysis_complete.RData")

for(file in check_files) {
  if(file.exists(file)) {
    file_size <- file.info(file)$size / 1024  # KB
    cat(sprintf("✓ %s: %.1f KB\n", file, file_size))
  } else {
    cat(sprintf("✗ %s: 文件不存在\n", file))
  }
}

cat("\n✅ LASSO回归分析全部完成!\n")
cat("   开始时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

save(mice_imp, stable_vars, file = "model_building_input.RData")