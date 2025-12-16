---
title: 'limma_py: A Python Implementation of limma for Differential Expression Analysis'
tags:
  - Python
  - limma
  - proteomics
  - differential expression
authors:
  - name: Xian Zhang
    affiliation: "1"
  - name: Hongchao Ji
    orcid: 0000-0002-7364-0741
    corresponding: true
    affiliation: "1"
affiliations:
  - name: Shenzhen Branch, Guangdong Laboratory for Lingnan Modern Agriculture, Genome Analysis Laboratory, Ministry of Agriculture and Rural Affairs, Agricultural Genomics Institute at Shenzhen, Chinese Academy of Agricultural Sciences, Shenzhen, 518120, China PR.
    index: 1
date: 15 December 2025
bibliography: paper.bib
repository: https://github.com/hcji/limma_py
---

# Summary

Python has increasingly replaced traditional R-based environments in both academic and industrial research, owing to its simplicity, readability, and the maturity of its scientific ecosystem [@cock_biopython_2009]. `limma_py` is a Python-native implementation of the core statistical framework of **limma**, a widely used method for differential expression and differential abundance analysis based on linear models with empirical Bayes variance moderation [@ritchie_limma_2015]. Originally developed for microarray-based transcriptomics, limma has proven broadly applicable across high-dimensional omics data, including RNA-seq, quantitative proteomics, metabolomics, and protein stability profiling.

`limma_py` reproduces results consistent with the original R implementation of limma while enabling fully Python-based analytical workflows. The package accepts standard Python data structures such as `pandas.DataFrame` and `numpy.ndarray`, and returns results in `DataFrame` format, facilitating seamless integration with downstream data analysis, visualization, and machine learning pipelines.

# Statement of Need

Quantitative proteomics experiments routinely rely on statistical comparisons across biological or treatment conditions to identify proteins that are differentially abundant or exhibit altered stability. Robust analysis of such data requires statistical models that can accommodate biological variability, small sample sizes, and the correction for multiple hypothesis testing. The linear modeling framework implemented in **limma** has become a reference method for these tasks [@ball_isothermal_2020], owing to its flexibility and the use of empirical Bayes variance moderation to improve statistical power and control false discovery rates.

Although originally developed for transcriptomic data, limma has demonstrated broad applicability across multiple omics modalities, including quantitative proteomics and thermal stability profiling assays. In particular, methods such as Thermal Proteome Profiling (TPP) and the Isothermal Shift Assay (iTSA) generate high-dimensional protein abundance measurements across experimental conditions, for which limma-style linear modeling provides a well-established and interpretable statistical solution.

At the same time, proteomics data analysis workflows have increasingly transitioned toward Python-based environments. Python has become a central platform for mass spectrometry data processing, machine learning, and integrative data analysis, supported by mature libraries such as `numpy`, `pandas`, `scikit-learn`, and proteomics-specific tools including `pyteomics` and `alphapept`. However, the original R-based implementation of limma presents practical barriers to integration within these Python-native workflows. While R–Python interfaces offer partial solutions, they introduce additional dependencies, data transfer overhead, and reduced transparency in parameter control.

`limma_py` addresses this gap by providing a Python-native reimplementation of the core limma methodology. By reproducing the statistical behavior of R limma with numerical equivalence while enabling direct interoperability with Python’s scientific ecosystem, `limma_py` allows proteomics researchers to perform differential abundance and stability analyses within a fully Python-based workflow. This supports reproducible, modular, and scalable analysis pipelines that align with current computational practices in quantitative proteomics.

# Usage Example

Below is a minimal example demonstrating how to perform differential expression analysis using `limma_py`. The workflow follows the standard limma procedure: defining an experimental design, fitting a linear model, specifying contrasts, applying empirical Bayes moderation, and extracting ranked results.

```python
import pandas as pd
import numpy as np
import limma_py

# Read CSV file, header=0 indicates the first row contains column names
data = pd.read_csv("data/Harmine_iTSA.csv", header=0)

# Extract gene names (first column)
gene_names = data.iloc[:, 0].values

# Extract expression matrix (all columns starting from the second)
expr_data = data.iloc[:, 1:]

# Define experimental groups: first 5 samples are control group (V_group), last 5 are treatment group (D_group)
group = np.array(["V_group"] * 5 + ["D_group"] * 5)

# Create design matrix using one-hot encoding
design_df = pd.get_dummies(group, drop_first=False)[["V_group", "D_group"]]

# Ensure design matrix is of integer type
design = design_df.astype(int)

# Create a copy of the expression matrix and set gene names as row indices
expr_matrix = expr_data.copy()
expr_matrix.index = gene_names

# Perform linear model fitting using limma
fit_python = limma_py.lmFit(expr_matrix, design)

# Set up contrast matrix: compare differences between treatment group (D_group) and control group (V_group)
contrasts = limma_py.make_contrasts('D_group - V_group', levels=design)

# Perform contrast analysis on the fitted results
fit_python = limma_py.contrasts_fit(fit_python, contrasts)

# Moderate standard errors using empirical Bayes method
eb_python = limma_py.eBayes(fit_python)

# Extract differential expression analysis result table
res = limma_py.toptable(eb_python)
```

# Validation on protein thermal stability datasets

To validate the correctness and practical applicability of `limma_py`, we compared its results with those obtained using the original R implementation of limma on four protein thermal stability datasets derived from the ProSAP study. These datasets include comparisons between control samples treated with dimethyl sulfoxide (DMSO) and drug-treated samples exposed to either 5-fluorouracil (5FU) or methotrexate (MTX), measured in two experimental contexts: intact cells and cell lysates [@ji_prosap_2022].

For each dataset, we evaluated differential protein abundance using identical experimental designs and contrast specifications, focusing on two aspects: statistical consistency and computational behavior. As shown in Figure 1, `limma_py` reproduces the results of R limma with numerical equivalence. Log fold change estimates and p-values are identical across all tested datasets, with all data points aligning precisely along the diagonal. This demonstrates that `limma_py` faithfully reproduces the statistical behavior of the original limma framework when applied to protein thermal stability data.

![Statistical consistency between `limma_py` and R limma on protein thermal stability datasets.(A) Comparison of log fold change estimates. (B) Comparison of p-values. Each point represents a protein, with values from `limma_py` plotted against those from R limma.](fig_1.png)

We further assessed computational performance by comparing the execution time of individual analysis steps. Overall runtimes of `limma_py` were comparable to those of R limma, with modest differences observed across specific functions. In particular, the linear model fitting step (`lmFit`) exhibited faster execution in Python, while empirical Bayes moderation (`eBayes`) showed similar performance between implementations. The result extraction step (`topTable`) was comparatively slower in Python. These results indicate that `limma_py` provides computational performance suitable for practical use in proteomics workflows while maintaining full statistical equivalence with the R implementation (Figure 2).

![Computational performance comparison between `limma_py` and R limma. (A) Execution time by analysis step. (B) Overall runtime comparison across datasets.](fig_2.png)

# Code availability

We have provided limma-py as a pip-installable package, available on PyPI at https://pypi.org/project/limma-py. The source code, including documentation and example notebooks, is openly available on GitHub at https://github.com/hcji/limma_py.

# Data availability

The dataset of benchmarking datasets can be found at https://github.com/hcji/ProSAP/tree/master/data.

# Acknowledgements

This work was supported by the Agricultural Science and Technology Innovation Program (CAAS-ZDRW202503)

# References
