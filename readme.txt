# 过敏性鼻炎神经阻断术疗效预测模型 - 分析代码

## 文件说明
- `01_data_prep_lasso.R`：数据预处理、多重插补、LASSO变量筛选
- `02_model_building.R`：多因素Logistic回归建模、模型验证
- `03_subgroup_analysis.R`：亚组分析、交互作用检验
- `shiny_app/app.R`：在线预测计算器源码
- `example_data.csv`：模拟示例数据（用于演示分析流程）
- `sessionInfo.txt`：R包版本信息

## 运行顺序
1. 先运行 `01_data_prep_lasso.R`
2. 再运行 `02_model_building.R`
3. 最后运行 `03_subgroup_analysis.R`

## 注意事项
- 原始数据因伦理要求无法公开，`example_data.csv` 为随机生成的模拟数据
- 本代码仅供学术验证，不得用于商业用途