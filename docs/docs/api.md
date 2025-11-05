# API参考

## 核心函数

### `lmFit(object, design=None, **kwargs)`
线性模型拟合函数，是limma分析流程的第一步。

**参数：**

- `object`: 基因表达数据矩阵，形状为(基因数, 样本数)，可以是numpy数组、pandas DataFrame

- `design`: 设计矩阵，形状为(样本数, 系数)，如果为None则创建仅截距模型

- `ndups`: 每个基因的重复数（用于同一基因由多个探针表示的情况）

- `spacing`: 重复之间的间距

- `weights`: 观测权重，可以是向量或与数据维度匹配的矩阵

- `method`: 拟合方法，目前仅支持"ls"（最小二乘法）


**返回：**
包含线性模型拟合结果的字典：
- `coefficients`: 系数矩阵 (基因数 × 系数)

- `stdev_unscaled`: 未缩放的标准误差矩阵

- `sigma`: 残差标准差

- `df_residual`: 残差自由度

- `cov_coefficients`: 系数协方差矩阵

- `rank`: 设计矩阵的秩

- `Amean`: 跨样本的平均表达水平

- `design`: 使用的设计矩阵


### `make_contrasts(*args, contrasts=None, levels)`
生成对比矩阵，支持直接输入pandas DataFrame作为levels参数。

**参数：**

- `*args`:用于传递对比表达式字符串（例如 "groupD - groupC"），与 contrasts 参数功能相同，两者可选择其一或同时使用
示例：make_contrasts("A - B", "C - (A+B)/2", levels=levels)

- `contrasts`: 用于传递对比表达式字符串（例如 "groupD - groupC"），与 contrasts 参数功能相同，两者可选择其一或同时使用
示例：make_contrasts("A - B", "C - (A+B)/2", levels=levels)

- `levels`: 设计矩阵（可直接输入pd.DataFrame）或列名列表

**返回：**
- `对比矩阵`: (np.ndarray)形状为 (n_levels, n_contrasts)，其中：
n_levels 是水平数量（来自 levels 参数）
n_contrasts 是对比表达式数量（来自 *args 和 contrasts 参数）

### `contrasts_fit(fit, contrasts=None, coefficients=None)`
从线性模型拟合结果中提取指定组对比的结果。

**参数：**

- `fit`: 线性模型拟合结果字典，必须包含以下键：
  
coefficients: 基因×系数表达系数矩阵
  
stdev_unscaled: 基因×系数未缩放标准差矩阵
  
cov_coefficients: 系数协方差矩阵（可选）

- `contrasts`: 对比矩阵（行数 = fit中的系数数量，列数 = 对比数量）

- `coefficients`: 要保留的系数列索引/名称（与contrasts互斥）

**返回：**
更新后的拟合结果字典，核心字段根据对比矩阵调整。

### `eBayes(fit, proportion=0.01, trend=False, **kwargs)`

线性模型拟合的经验贝叶斯调节标准误差。

**参数：**
- `fit`: 来自lmFit的线性模型拟合结果字典

- `proportion`: 预期差异表达基因的先验比例，默认0.01

- `trend`: 是否考虑方差随表达水平的趋势

- `winsor_tail_p`: 稳健估计中的Winsorization尾概率


**返回：**
更新后的拟合字典，包含经验贝叶斯结果：
- `df_prior`: 先验自由度

- `s2_prior`: 先验方差

- `t`: 调节后的t统计量

- `p_value`: 调节后的p值

- `lods`: 差异表达的对数概率

- `F`: F统计量（如果设计允许）

- `F_p_value`: F检验p值

### `toptable(fit, coef=1, number=10, adjust_method="BH", **kwargs)`

从线性模型拟合中提取排名靠前的基因表。

**参数：**
- `fit`: 通常来自eBayes的线性模型拟合结果字典

- `coef`: 要显示的系数，可以是列索引或索引列表

- `number`: 要显示的顶级基因数量，默认10

- `adjust_method`: 多重检验校正方法："BH"（Benjamini-Hochberg）、"BY"、"bonferroni"、"holm"、"none"

- `sort_by`: 结果排序依据："B"（对数概率）、"logFC"、"AveExpr"、"P"（p值）、"t"、"none"

- `p_value`: p值过滤 cutoff，默认1.0（无过滤）

- `lfc`: 对数倍变化cutoff（log2尺度）

- `confint`: 是否计算置信区间


**返回：**
包含顶级基因的pandas DataFrame，列包括：
- `logFC`: 对数倍变化

- `AveExpr`: 平均表达水平

- `t`: 调节后的t统计量

- `P.Value`: p值

- `adj.P.Val`: 调整后的p值（FDR）

- `B`: 差异表达的对数概率

- `CI.L`, `CI.R`: 置信区间上下限（如果confint=True）
