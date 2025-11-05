## 数据格式要求
limma-py要求输入数据为CSV格式，具体结构如下：
```csv
Gene,Sample1,Sample2,Sample3,Sample4,Sample5,Sample6,Sample7,Sample8,Sample9,Sample10
Gene_1,8.45,7.89,8.12,15.67,14.89,16.23,8.33,7.95,8.21,8.09
Gene_2,12.34,11.89,12.56,11.45,10.98,12.11,25.67,24.89,26.01,25.34
Gene_3,5.67,6.01,5.89,18.90,19.45,18.23,6.12,5.78,6.34,5.95
```

### 格式说明：

- 第一列: 基因/蛋白名称（任意标识符）

- 后续列: 数值型表达矩阵，每列代表一个样本

- 必须为CSV格式，其他格式需要先转换

## 实际示例
```python
import pandas as pd
import numpy as np
import limma_py

# 读取CSV文件，header=0表示第一行是列名
data = pd.read_csv("data/Harmine_iTSA.csv", header=0)

# 提取基因名称（第一列）
gene_names = data.iloc[:, 0].values

# 提取表达矩阵（从第二列开始的所有列）
expr_data = data.iloc[:, 1:]

# 定义实验组别：前5个样本为对照组(V_group)，后5个样本为处理组(D_group)
group = np.array(["V_group"] * 5 + ["D_group"] * 5)

# 使用one-hot编码创建设计矩阵
design_df = pd.get_dummies(group, drop_first=False)[["V_group", "D_group"]]

# 确保设计矩阵为整数类型
design = design_df.astype(int)

# 创建表达矩阵副本，并设置基因名为行索引
expr_matrix = expr_data.copy()
expr_matrix.index = gene_names

# 使用limma进行线性模型拟合
fit_python = limma_py.lmFit(expr_matrix, design)

# 设置对比矩阵：比较处理组(D_group)与对照组(V_group)的差异
contrasts = limma_py.make_contrasts('D_group - V_group', levels=design)

# 对拟合结果进行对比分析
fit_python = limma_py.contrasts_fit(fit_python, contrasts)

# 使用经验贝叶斯方法调节标准误差
eb_python = limma_py.eBayes(fit_python)

# 提取差异表达分析结果表格
res = limma_py.toptable(eb_python)
```
## 常见问题
❌ 错误提示："数据维度不匹配"
- 检查设计矩阵的行数是否与表达数据的列数（样本数）一致

- 检查组别向量的长度是否与样本数一致

✅ 数据预处理检查清单
- CSV文件第一列为基因名

- 第一行为样本名

- 所有表达值为数值型

- 无缺失值或已处理缺失值

- 样本顺序与组别定义一致