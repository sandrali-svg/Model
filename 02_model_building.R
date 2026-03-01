# ==================== 初始化设置 ====================
rm(list = ls())
gc()

if(!file.exists("model_building_input.RData")) {
  stop("请先运行 01_data_prep_lasso.R 生成所需数据文件")
}
load("model_building_input.RData")

# 加载必要的包
required_packages <- c("mice", "rms", "boot", "pROC", "caret", "ggplot2", 
                       "openxlsx", "gridExtra", "ggpubr", "reshape2", "dplyr",
                       "rmda", "Hmisc", "ResourceSelection", "parallel")

for(pkg in required_packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

set.seed(123456)

# ==================== 数据准备 ====================
cat("=== 数据准备 ===\n")

# 加载LASSO筛选结果（应包含 mice_imp 对象）
load("model_building_input_filtered.RData")

# 检查是否包含多重插补对象
if(!exists("mice_imp")) stop("未找到 mice_imp 对象")

# 定义数据转换函数（用于构建 imp_list）
transform_data <- function(d) {
  d$疗效_num <- ifelse(d$疗效 == "显效", 1, 0)
  d$LMR <- d$`淋巴细胞与单核细胞的比率`
  d$phadiatop试验分级 <- factor(d$`吸入变应原筛查分级`, levels = 0:6)
  d$是否伴OSA <- factor(d$是否伴OSA, levels = c(0, 1), labels = c("否", "是"))
  return(d)
}

# 将转换应用到所有插补数据集（生成完整数据列表）
imp_list <- complete(mice_imp, "long")  # 长格式，包含 .imp 和 .id
imp_list <- split(imp_list, imp_list$.imp)
imp_list <- lapply(imp_list, transform_data)

# 样本量信息（以第一个数据集为例）
sample_data <- imp_list[[1]]
cat("总样本量:", nrow(sample_data), "\n")
cat("显效(事件):", sum(sample_data$疗效_num == 1), "例\n")
cat("有效(非事件):", sum(sample_data$疗效_num == 0), "例\n")
cat("事件发生率:", round(mean(sample_data$疗效_num) * 100, 1), "%\n")

# ==================== 多因素逻辑回归（合并100个插补数据集） ====================
cat("\n=== 多因素逻辑回归（合并100个插补数据集） ===\n")

# 定义公式（使用转换后的变量名）
full_formula <- 疗效_num ~ 术前VAS评分 + 病程 + LMR + 空腹葡萄糖 + 是否伴OSA + phadiatop试验分级

# 对每个插补数据集拟合模型（在 with 内部直接转换变量）
models <- with(mice_imp, {
  疗效_num <- ifelse(疗效 == "显效", 1, 0)
  LMR <- `淋巴细胞与单核细胞的比率`
  phadiatop试验分级 <- factor(`吸入变应原筛查分级`, levels = 0:6)
  是否伴OSA <- factor(是否伴OSA, levels = c(0, 1), labels = c("否", "是"))
  glm(疗效_num ~ 术前VAS评分 + 病程 + LMR + 空腹葡萄糖 + 是否伴OSA + phadiatop试验分级,
      family = binomial)
})

# 提取每个插补数据集的模型对象列表
model_list <- models$analyses

# 合并结果
pooled <- pool(models)
pooled_summary <- summary(pooled, conf.int = TRUE, exponentiate = TRUE)

cat("\n合并后的多因素逻辑回归完成\n")

# ==================== 计算个体平均预测概率 ====================
cat("\n=== 计算个体平均预测概率 ===\n")

# 利用每个插补数据集对应的模型，预测该数据集本身的概率，然后平均
pred_probs_list <- lapply(1:length(imp_list), function(i) {
  predict(model_list[[i]], newdata = imp_list[[i]], type = "response")
})
pred_mat <- do.call(cbind, pred_probs_list)
avg_pred_prob <- rowMeans(pred_mat)

# ==================== ROC 曲线（基于平均预测概率） ====================
cat("\n=== 绘制 ROC 曲线 ===\n")

roc_obj <- roc(sample_data$疗效_num, avg_pred_prob)
auc_value <- auc(roc_obj)
auc_ci <- ci.auc(roc_obj)

# 核心新增：计算Youden指数最佳截断值
# best.method = "youden" 指定用Youden指数找最优截断点
best_cutoff <- coords(roc_obj, x = "best", best.method = "youden", 
                      ret = c("threshold", "sensitivity", "specificity"))

# 格式化结果：截断值保留3位小数，敏感度/特异度转换为百分比并保留1位小数
cutoff_value <- round(best_cutoff$threshold, 3)  # 最佳截断值（预测概率）
sensitivity_pct <- round(best_cutoff$sensitivity * 100, 1)  # 敏感度（百分比）
specificity_pct <- round(best_cutoff$specificity * 100, 1)  # 特异度（百分比）

# 输出你需要的结果语句
cat("\n=== Youden指数最优截断值分析结果 ===\n")
result_sentence <- paste0("根据Youden指数确定的最佳截断值为预测概率", cutoff_value, 
                          "，此时模型的敏感度为", sensitivity_pct, "%，特异度为", 
                          specificity_pct, "%，意味着该模型能正确识别出", sensitivity_pct, 
                          "%的显效患者，且误判率较低。")
cat(result_sentence, "\n")

# 绘制并保存ROC曲线（保留你原有绘图逻辑）
png("ROC_curve_术后疗效预测模型.png", width = 8, height = 8, units = "in", res = 300)
par(mar = c(4,4,1,1))
plot(roc_obj, col = "#2E8B57", lwd = 3, print.auc = TRUE, print.auc.x = 0.8, print.auc.y = 0.1,
     print.auc.cex = 1.2, grid = TRUE, grid.col = "gray90", legacy.axes = TRUE,
     xlab = "1-特异度", ylab = "敏感度", main = "")
legend("bottomright", legend = paste0("AUC = ", round(auc_value,3), " (95%CI: ", 
                                      round(auc_ci[1],3), "-", round(auc_ci[3],3), ")"),
       col = "#2E8B57", lwd = 3, bty = "n", cex = 1.1)
dev.off()
cat("ROC曲线已保存\n")

# ==================== 生成最终模型回归结果表（表6） ====================
cat("\n=== 生成最终模型回归结果表（表6） ===\n")

# 从pooled_summary提取数据，排除截距项
pooled_res <- pooled_summary[pooled_summary$term != "(Intercept)", ]

# 计算β系数（log OR）并保留四位小数
pooled_res$beta <- log(pooled_res$estimate)

# 构建表格所需列
table6 <- data.frame(
  变量 = pooled_res$term,
  β系数 = round(pooled_res$beta, 4),
  标准误 = round(pooled_res$std.error, 4),
  P值 = round(pooled_res$p.value, 4),
  OR值 = round(pooled_res$estimate, 4),
  OR_95CI_lower = round(pooled_res$`2.5 %`, 4),
  OR_95CI_upper = round(pooled_res$`97.5 %`, 4)
)

# 查看表格
print(table6)

# 保存为CSV和Excel
write.csv(table6, "最终模型回归结果_表6.csv", row.names = FALSE, fileEncoding = "UTF-8")
library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, "表6 多因素回归结果")
writeData(wb, "表6 多因素回归结果", table6)
saveWorkbook(wb, "最终模型回归结果_表6.xlsx", overwrite = TRUE)

cat("表6已保存为 CSV 和 Excel 格式。\n")

# ==================== Bootstrap 内部验证（强制包含所有因子水平） ====================
cat("\n=== Bootstrap内部验证（1000次，强制包含所有因子水平） ===\n")

# 获取所有插补数据集中因子的完整水平（并集）
all_phadiatop_levels <- unique(unlist(lapply(imp_list, function(df) levels(df$phadiatop试验分级))))
all_osa_levels <- unique(unlist(lapply(imp_list, function(df) levels(df$是否伴OSA))))

n_boot <- 1000
boot_auc <- numeric(n_boot)
pb <- txtProgressBar(min = 0, max = n_boot, style = 3)

for(i in 1:n_boot) {
  # 随机选取一个插补数据集
  imp_idx <- sample(1:length(imp_list), 1)
  boot_data <- imp_list[[imp_idx]]
  
  # 确保当前数据集的因子水平为完整水平
  boot_data$phadiatop试验分级 <- factor(boot_data$phadiatop试验分级, levels = all_phadiatop_levels)
  boot_data$是否伴OSA <- factor(boot_data$是否伴OSA, levels = all_osa_levels)
  
  # 常规Bootstrap抽样（有放回）
  boot_idx <- sample(1:nrow(boot_data), nrow(boot_data), replace = TRUE)
  boot_sample <- boot_data[boot_idx, ]
  
  # 检查每个因子水平的频数，若缺失则从原始数据中补充
  missing_levels <- setdiff(all_phadiatop_levels, levels(droplevels(boot_sample$phadiatop试验分级)))
  if(length(missing_levels) > 0) {
    for(lvl in missing_levels) {
      # 从原始数据中随机抽取一条该水平的观测
      add_row <- boot_data[boot_data$phadiatop试验分级 == lvl, ][1, ]  # 取第一条
      if(!is.na(add_row[1,1])) {
        boot_sample <- rbind(boot_sample, add_row)
      }
    }
  }
  
  # 再次确保因子水平
  boot_sample$phadiatop试验分级 <- factor(boot_sample$phadiatop试验分级, levels = all_phadiatop_levels)
  boot_sample$是否伴OSA <- factor(boot_sample$是否伴OSA, levels = all_osa_levels)
  
  # 在Bootstrap样本上拟合模型
  boot_fit <- tryCatch({
    glm(full_formula, data = boot_sample, family = binomial)
  }, error = function(e) NULL)
  
  if(is.null(boot_fit) || any(is.na(coef(boot_fit)))) {
    boot_auc[i] <- NA
    setTxtProgressBar(pb, i)
    next
  }
  
  # 用该模型预测所有插补数据集的概率，并取平均
  pred_boot <- sapply(imp_list, function(df) {
    # 预测前，将每个数据集的因子水平对齐到完整水平
    df$phadiatop试验分级 <- factor(df$phadiatop试验分级, levels = all_phadiatop_levels)
    df$是否伴OSA <- factor(df$是否伴OSA, levels = all_osa_levels)
    predict(boot_fit, newdata = df, type = "response")
  })
  pred_boot_avg <- rowMeans(pred_boot)
  
  # 计算 AUC
  boot_auc[i] <- as.numeric(roc(sample_data$疗效_num, pred_boot_avg)$auc)
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 移除 NA 值
valid_boot_auc <- boot_auc[!is.na(boot_auc)]
valid_count <- length(valid_boot_auc)

if(valid_count < 500) {
  cat("\n警告：有效Bootstrap样本不足500次（实际", valid_count, "次），结果可能不稳定\n")
}

# 计算乐观校正
original_auc <- auc_value
optimism <- mean(valid_boot_auc) - original_auc
corrected_auc <- original_auc - optimism

cat("\nBootstrap验证结果:\n")
cat("  有效Bootstrap样本数:", valid_count, "/", n_boot, "\n")
cat("  原始C-statistic (AUC):", round(original_auc, 3), "\n")
cat("  Bootstrap平均AUC:", round(mean(valid_boot_auc), 3), "\n")
cat("  乐观估计:", round(optimism, 3), "\n")
cat("  校正后C-statistic:", round(corrected_auc, 3), "\n")
cat("  95%置信区间 (Bootstrap百分位数):", round(quantile(valid_boot_auc, c(0.025, 0.975)), 3), "\n")

# ==================== 决策曲线（基于平均预测概率） ====================
cat("\n=== 绘制决策曲线 ===\n")

# 准备数据框（包含结局和平均预测概率）
dc_data <- data.frame(疗效_num = sample_data$疗效_num, pred = avg_pred_prob)

tryCatch({
  dc_obj <- decision_curve(疗效_num ~ pred, data = dc_data,
                           family = binomial(link = "logit"),
                           thresholds = seq(0, 1, by = 0.01),
                           confidence.intervals = 0.95,
                           study.design = "cohort")
  
  png("decision_curve_术后疗效预测模型.png", width = 8, height = 8, units = "in", res = 300)
  par(mar = c(4,4,1,1))
  plot_decision_curve(dc_obj, curve.names = "预测模型", cost.benefit.axis = FALSE,
                      col = "#8E44AD", confidence.intervals = FALSE, standardize = FALSE,
                      legend.position = "topright", lwd = 3,
                      xlab = "阈值概率", ylab = "净收益", main = "")
  grid(col = "gray90", lty = "dotted")
  dev.off()
  cat("决策曲线已保存\n")
}, error = function(e) {
  cat("决策曲线生成失败:", e$message, "\n")
})

# ==================== 列线图构建（使用 fit.mult.impute 合并模型，实时转换数据） ====================
cat("\n=== 构建列线图（基于 fit.mult.impute 合并模型） ===\n")

# 从第一个插补数据集生成模板数据（转换后，用于 datadist）
template_data <- transform_data(complete(mice_imp, 1))
dd <- datadist(template_data)
options(datadist = "dd")

# 定义拟合函数：先转换数据，再拟合 lrm
fit_lrm <- function(formula, data) {
  data <- transform_data(data)   # 转换变量
  lrm(formula, data = data, x = TRUE, y = TRUE)
}

# 使用 fit.mult.impute 合并模型（直接传入 mice_imp 对象）
fmi <- fit.mult.impute(full_formula, fitter = fit_lrm, xtrans = mice_imp, data = template_data)

# 绘制列线图
nom <- nomogram(fmi, 
                fun = function(x) plogis(x), 
                fun.at = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95),
                funlabel = "预测概率",
                lp = FALSE,
                conf.int = FALSE)

png("nomogram_术后疗效预测模型.png", width = 10, height = 8, units = "in", res = 300)
par(mar = c(4,4,1,1))
plot(nom, cex.axis = 0.7, cex.var = 0.8, lmgp = 0.3, xfrac = 0.15)
dev.off()
cat("列线图已保存（基于 fit.mult.impute 合并模型）\n")

# ==================== 校准曲线（基于合并模型，Bootstrap 1000次） ====================
cat("\n=== 绘制校准曲线（Bootstrap 1000次） ===\n")

# 对合并后的模型进行校准验证
cal <- calibrate(fmi, method = "boot", B = 1000)

png("calibration_curve_术后疗效预测模型.png", width = 8, height = 8, units = "in", res = 300)
par(mar = c(4,4,1,1))
plot(cal,
     xlim = c(0, 1),
     ylim = c(0, 1),
     xlab = "预测概率",
     ylab = "实际概率",
     subtitles = FALSE,
     main = "")
abline(0, 1, col = "red", lty = 2, lwd = 2)
dev.off()
cat("校准曲线已保存（Bootstrap 1000次）\n")

# ==================== 最终总结 ====================
cat("\n========== 分析完成总结 ==========\n")
cat("\n1. 样本特征:\n")
cat("   总样本量:", nrow(sample_data), "例\n")
cat("   显效(事件):", sum(sample_data$疗效_num == 1), "例\n")
cat("   有效(非事件):", sum(sample_data$疗效_num == 0), "例\n")
cat("   事件发生率:", round(mean(sample_data$疗效_num) * 100, 1), "%\n")

cat("\n2. 模型构建:\n")
cat("   纳入变量: 6个（基于LASSO稳定性筛选）\n")
cat("   建模方法: 多因素logistic回归（合并100个插补数据集）\n")

cat("\n3. 模型性能（基于平均预测概率）:\n")
cat("   C-index (AUC):", round(auc_value, 3), 
    " (", round(auc_ci[1], 3), "-", round(auc_ci[3], 3), ")\n")
cat("   校正后C-index (Bootstrap 1000次):", round(corrected_auc, 3), "\n")

cat("\n4. 生成文件清单:\n")
cat("   可视化图形:\n")
cat("     - ROC_curve_术后疗效预测模型.png\n")
cat("     - calibration_curve_术后疗效预测模型.png\n")
cat("     - decision_curve_术后疗效预测模型.png\n")
cat("     - nomogram_术后疗效预测模型.png\n")

cat("\n✅ 术后疗效预测模型分析（整合100个插补数据集）全部完成！\n")
cat("   当前时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

# ==================== 比较三个模型：模型A、B、C（基于平均预测概率） ====================
cat("\n=== 比较模型A、B、C性能（基于平均预测概率） ===\n")

# 定义三个模型的公式
formula_A <- 疗效_num ~ 病程 + 空腹葡萄糖 + 是否伴OSA + 术前VAS评分
formula_B <- 疗效_num ~ 病程 + 空腹葡萄糖 + 是否伴OSA + 术前VAS评分 + phadiatop试验分级
formula_C <- 疗效_num ~ 术前VAS评分 + 病程 + LMR + 空腹葡萄糖 + 是否伴OSA + phadiatop试验分级

n_imp <- length(imp_list)
n_samples <- nrow(imp_list[[1]])

# 存储每个模型在每个数据集上的预测概率矩阵（行=个体，列=插补数据集）
pred_probs_A <- matrix(NA, nrow = n_samples, ncol = n_imp)
pred_probs_B <- matrix(NA, nrow = n_samples, ncol = n_imp)
pred_probs_C <- matrix(NA, nrow = n_samples, ncol = n_imp)

# 存储每个模型的AIC（每个数据集一个值）
AIC_A <- numeric(n_imp)
AIC_B <- numeric(n_imp)
AIC_C <- numeric(n_imp)

# 真实结局（所有插补数据集结局相同，用第一个即可）
y_true <- imp_list[[1]]$疗效_num

cat("正在对100个插补数据集分别拟合三个模型并收集预测概率...\n")
pb <- txtProgressBar(min = 1, max = n_imp, style = 3)

for(k in 1:n_imp) {
  dat <- imp_list[[k]]
  
  # 模型A
  fitA <- glm(formula_A, data = dat, family = binomial)
  pred_probs_A[, k] <- predict(fitA, newdata = dat, type = "response")
  AIC_A[k] <- AIC(fitA)
  
  # 模型B
  fitB <- glm(formula_B, data = dat, family = binomial)
  pred_probs_B[, k] <- predict(fitB, newdata = dat, type = "response")
  AIC_B[k] <- AIC(fitB)
  
  # 模型C
  fitC <- glm(formula_C, data = dat, family = binomial)
  pred_probs_C[, k] <- predict(fitC, newdata = dat, type = "response")
  AIC_C[k] <- AIC(fitC)
  
  setTxtProgressBar(pb, k)
}
close(pb)

# 计算每个个体的平均预测概率
avg_pred_A <- rowMeans(pred_probs_A)
avg_pred_B <- rowMeans(pred_probs_B)
avg_pred_C <- rowMeans(pred_probs_C)

# 计算每个模型的AUC（基于平均预测概率）及其95% DeLong置信区间
library(pROC)
roc_A <- roc(y_true, avg_pred_A)
auc_A <- auc(roc_A)
auc_ci_A <- ci.auc(roc_A)

roc_B <- roc(y_true, avg_pred_B)
auc_B <- auc(roc_B)
auc_ci_B <- ci.auc(roc_B)

roc_C <- roc(y_true, avg_pred_C)
auc_C <- auc(roc_C)
auc_ci_C <- ci.auc(roc_C)

# 计算敏感性和特异性（阈值0.5）
threshold <- 0.5
calc_sens_spec <- function(pred, true) {
  pred_class <- ifelse(pred > threshold, 1, 0)
  tp <- sum(pred_class == 1 & true == 1)
  tn <- sum(pred_class == 0 & true == 0)
  fp <- sum(pred_class == 1 & true == 0)
  fn <- sum(pred_class == 0 & true == 1)
  sens <- tp / (tp + fn)
  spec <- tn / (tn + fp)
  return(c(sens, spec))
}

sens_spec_A <- calc_sens_spec(avg_pred_A, y_true)
sens_spec_B <- calc_sens_spec(avg_pred_B, y_true)
sens_spec_C <- calc_sens_spec(avg_pred_C, y_true)

# 对于AIC，计算100个数据集的均值及95%百分位数置信区间
calc_aic_ci <- function(x) {
  c(mean = mean(x), lower = quantile(x, 0.025), upper = quantile(x, 0.975))
}
aic_summary_A <- calc_aic_ci(AIC_A)
aic_summary_B <- calc_aic_ci(AIC_B)
aic_summary_C <- calc_aic_ci(AIC_C)

# 构建表5
table5 <- data.frame(
  模型 = c("模型A", "模型B", "模型C"),
  变量数 = c(4, 5, 6),
  AUC = c(
    sprintf("%.3f (%.3f-%.3f)", auc_A, auc_ci_A[1], auc_ci_A[3]),
    sprintf("%.3f (%.3f-%.3f)", auc_B, auc_ci_B[1], auc_ci_B[3]),
    sprintf("%.3f (%.3f-%.3f)", auc_C, auc_ci_C[1], auc_ci_C[3])
  ),
  AIC = c(
    sprintf("%.1f (%.1f-%.1f)", aic_summary_A["mean"], aic_summary_A["lower"], aic_summary_A["upper"]),
    sprintf("%.1f (%.1f-%.1f)", aic_summary_B["mean"], aic_summary_B["lower"], aic_summary_B["upper"]),
    sprintf("%.1f (%.1f-%.1f)", aic_summary_C["mean"], aic_summary_C["lower"], aic_summary_C["upper"])
  ),
  敏感性 = sprintf("%.3f", c(sens_spec_A[1], sens_spec_B[1], sens_spec_C[1])),
  特异性 = sprintf("%.3f", c(sens_spec_A[2], sens_spec_B[2], sens_spec_C[2]))
)

# 设置列名
names(table5) <- c("模型", "变量数", "AUC (95% CI)", "AIC (95% CI)", "敏感性", "特异性")

# 打印表格
cat("\n========== 表5 模型性能比较 ==========\n")
print(table5, row.names = FALSE)

# 保存为CSV和Excel
write.csv(table5, "模型性能比较_表5.csv", row.names = FALSE, fileEncoding = "UTF-8")
library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, "表5 模型性能比较")
writeData(wb, "表5 模型性能比较", table5)
saveWorkbook(wb, "模型性能比较_表5.xlsx", overwrite = TRUE)
cat("\n表5已保存为 CSV 和 Excel 格式。\n")

# ==================== 生成最终模型回归结果表（表6） ====================
cat("\n=== 生成最终模型回归结果表（表6） ===\n")

# 从pooled_summary提取数据，排除截距项
pooled_res <- pooled_summary[pooled_summary$term != "(Intercept)", ]

# 计算β系数（log OR）并保留四位小数
pooled_res$beta <- log(pooled_res$estimate)

# 构建表格所需列
table6 <- data.frame(
  变量 = pooled_res$term,
  β系数 = round(pooled_res$beta, 4),
  标准误 = round(pooled_res$std.error, 4),
  P值 = round(pooled_res$p.value, 4),
  OR值 = round(pooled_res$estimate, 4),
  OR_95CI_lower = round(pooled_res$`2.5 %`, 4),
  OR_95CI_upper = round(pooled_res$`97.5 %`, 4)
)

# 查看表格
print(table6)

# 保存为CSV和Excel
write.csv(table6, "最终模型回归结果_表6.csv", row.names = FALSE, fileEncoding = "UTF-8")
library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, "表6 多因素回归结果")
writeData(wb, "表6 多因素回归结果", table6)
saveWorkbook(wb, "最终模型回归结果_表6.xlsx", overwrite = TRUE)

cat("表6已保存为 CSV 和 Excel 格式。\n")

# ==================== 检查二分类模型在各插补数据集中的方向一致性 ====================
cat("\n=== 检查二分类模型在各插补数据集中的方向一致性 ===\n")

# 从 model_list 中提取每个模型的系数符号
coef_sign_list <- lapply(model_list, function(fit) {
  coefs <- coef(fit)
  sign(coefs)  # 返回 -1, 0, 1
})

# 转换为矩阵，方便统计
sign_matrix <- do.call(rbind, coef_sign_list)

# 计算每个变量在所有数据集中系数为正的比例
pos_prop <- colMeans(sign_matrix == 1, na.rm = TRUE)
neg_prop <- colMeans(sign_matrix == -1, na.rm = TRUE)
zero_prop <- colMeans(sign_matrix == 0, na.rm = TRUE)

direction_check <- data.frame(
  变量 = colnames(sign_matrix),
  正方向比例 = round(pos_prop * 100, 1),
  负方向比例 = round(neg_prop * 100, 1),
  零比例 = round(zero_prop * 100, 1),
  结论 = ifelse(pos_prop >= 0.95, "绝大多数为正", 
              ifelse(neg_prop >= 0.95, "绝大多数为负", 
                     ifelse(pos_prop >= 0.8, "多数为正",
                            ifelse(neg_prop >= 0.8, "多数为负", "方向不稳定"))))
)

print(direction_check)

# 导出为Excel文件
if (!require(openxlsx)) install.packages("openxlsx")
library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, "方向一致性检查")
writeData(wb, "方向一致性检查", direction_check)
saveWorkbook(wb, "二分类方向一致性检查.xlsx", overwrite = TRUE)

cat("\n✅ 方向一致性检查结果已保存为 '二分类方向一致性检查.xlsx'\n")

# 方法二：从 pooled_summary 中提取系数（log OR）
if (exists("pooled_summary")) {
  # 提取变量名和 OR
  var_names <- pooled_summary$term
  OR_values <- pooled_summary$estimate
  
  # 计算 log OR
  logOR_values <- log(OR_values)
  
  # 打印结果
  cat("\n========== 模型系数（log OR）==========\n")
  for (i in seq_along(var_names)) {
    cat(sprintf("%-30s = %f\n", var_names[i], logOR_values[i]))
  }
} else {
  cat("错误：未找到 'pooled_summary' 对象，请先运行建模代码中的合并部分。\n")
}

