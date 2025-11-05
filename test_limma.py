import pandas as pd
import numpy as np
import limma_py

data = pd.read_csv("C:\\Users\\29675\\Documents\\Harmine_iTSA.csv", header=0)
gene_names = data.iloc[:, 0].values
expr_data = data.iloc[:, 1:]
group = np.array(["V_group"] * 5 + ["D_group"] * 5)

design_df = pd.get_dummies(group, drop_first=False)[["V_group", "D_group"]]

design = design_df.astype(int)
expr_matrix = expr_data.copy()
expr_matrix.index = gene_names
fit_python = limma_py.lmFit(expr_matrix, design)
contrasts = limma_py.make_contrasts('D_group - V_group',levels = design)
fit_python = limma_py.contrasts_fit(fit_python, contrasts)
eb_python = limma_py.eBayes(fit_python)
res = limma_py.toptable(eb_python,coef=[0])
