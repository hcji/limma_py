# limma_py

A Python implementation of the popular R limma package for differential expression analysis of gene expression data.

[![Python Version](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

## Overview

`limma_py` is a comprehensive Python port of the widely-used R limma (Linear Models for Microarray Data) package. It provides powerful statistical methods for analyzing gene expression data from microarray and RNA-seq experiments, with a focus on differential expression analysis.

### Key Features

- **Linear Model Fitting**: Fit linear models for each gene in expression data
- **Empirical Bayes Moderation**: Borrow information across genes to improve variance estimates
- **Contrast Analysis**: Handle complex experimental designs with multiple comparisons
- **Multiple Testing Correction**: Comprehensive p-value adjustment methods
- **Comprehensive Output**: Detailed statistical results with fold changes, p-values, and confidence intervals

## Installation

```bash
pip install limma_py