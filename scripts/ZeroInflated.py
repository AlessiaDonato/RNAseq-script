#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 12 15:20:06 2025

@author: donato.alessia
"""


#Identificare geni zero-inflated (ZI) — cioè geni che mostrano un numero anomalo 
#di zeri nelle matrici di espressione genica a singola cellula — in specifici sottotipi cellulari. 
#Lo fa utilizzando scVI-AutoZI, un modello probabilistico bayesiano che distingue 
#tra geni "realmente silenziati" e "falsi zeri" dovuti a dropout tecnici in scRNA-seq.
cd ..//..//beegfs/scratch/ric.boletta/donato.alessia/epithilial/scvl

#######Identification of zero-inflated genes

import tempfile

import numpy as np
import scanpy as sc
import scvi
import seaborn as sns
import torch
from scipy.stats import beta

import os
import anndata as ad
import scanpy as sc
import pandas as pd
import scanpy.external as sce
import schist as scs
import matplotlib as plt
plt.use('tkagg')

merged= ad.read_h5ad("/beegfs/scratch/ric.boletta/ric.boletta/tsc1/RCC_integrated/python/harmony_P33_P80_integrated_hvg_3000_no_outliers_n10.h5ad")


transforming_IC= merged[merged.obs["celltype"].isin(["KO_IC", "KO_transitional_regular",'transforming_1', 'transforming_2', 
                                                     'transforming_3','transforming_4', 'transforming_5'])]

transforming_IC = transforming_IC[transforming_IC.obs.pct_counts_mt <= 30, :]
transforming_IC.obs = transforming_IC.obs.iloc[:,:30]

transforming_IC
#View of AnnData object with n_obs × n_vars = 8593 × 27836

print(transforming_IC.obs["celltype"].value_counts())
#celltype
#transforming_2             1718
#KO_transitional_regular    1456
#transforming_4             1271
#transforming_5             1226
#transforming_1             1066
#KO_IC                      1014
#transforming_3              842

del transforming_IC.layers
transforming_IC.layers["counts"]= transforming_IC.X.copy()

transforming_IC.layers["counts"]= transforming_IC.X.copy()

scvi.data.poisson_gene_selection(
    transforming_IC,
    n_top_genes=3000,
    batch_key="batch",
    subset=True,
    layer="counts",
)

scvi.model.AUTOZI.setup_anndata(
    transforming_IC,
    labels_key="celltype",
    batch_key="batch",
    layer="counts",
)
model = scvi.model.AUTOZI(transforming_IC)
model.train(max_epochs=200, plan_kwargs={"lr": 1e-2})

outputs = model.get_alphas_betas()
alpha_posterior = outputs["alpha_posterior"]
beta_posterior = outputs["beta_posterior"]

# Threshold (or Kzinb/Knb+Kzinb in paper)
threshold = 0.5

# q(delta_g < 0.5) probabilities
zi_probs = beta.cdf(0.5, alpha_posterior, beta_posterior)

# ZI genes
is_zi_pred = zi_probs > threshold

print("Fraction of predicted ZI genes :", is_zi_pred.mean())
threshold = 0.5

# q(delta_g < 0.5) probabilities
zi_probs = beta.cdf(0.5, alpha_posterior, beta_posterior)
 
# ZI genes
is_zi_pred = zi_probs > threshold

print("Fraction of predicted ZI genes :", is_zi_pred.mean())
#Fraction of predicted ZI genes : 0.795
threshold = 0.75
# q(delta_g < 0.75) probabilities
is_zi_pred = zi_probs > threshold
print("Fraction of predicted ZI genes :", is_zi_pred.mean())
#Fraction of predicted ZI genes : 0.7516666666666667
threshold = 0.8
# q(delta_g < 0.8) probabilities
is_zi_pred = zi_probs > threshold
print("Fraction of predicted ZI genes :", is_zi_pred.mean())
#Fraction of predicted ZI genes : 0.7376666666666667


Analyze gene-cell-type-specific ZI

One may argue that zero-inflation should also be treated on the cell-type (or ‘label’) level, in addition to the gene level. AutoZI can be extended by assuming a random variable 
 for each gene 
 and cell type 
 which denotes the probability that gene 
 is not zero-inflated in cell-type 
. The analysis above can be extended to this new scale.

# Model definition
model_genelabel = scvi.model.AUTOZI(transforming_IC, dispersion="gene-label", zero_inflation="gene-label")

# Training
model_genelabel.train(max_epochs=200, plan_kwargs={"lr": 1e-2})

# Retrieve posterior distribution parameters
outputs_genelabel = model_genelabel.get_alphas_betas()
alpha_posterior_genelabel = outputs_genelabel["alpha_posterior"]
beta_posterior_genelabel = outputs_genelabel["beta_posterior"]

# q(delta_g < 0.5) probabilities
zi_probs_genelabel = beta.cdf(0.5, alpha_posterior_genelabel, beta_posterior_genelabel)

# ZI gene-cell-types
is_zi_pred_genelabel = zi_probs_genelabel > threshold

ct = transforming_IC.obs.celltype.astype("category")
codes = np.unique(ct.cat.codes)
cats = ct.cat.categories
for ind_cell_type, cell_type in zip(codes, cats, ):
    is_zi_pred_genelabel_here = is_zi_pred_genelabel[:, ind_cell_type]
    print(
        f"Fraction of predicted ZI genes for cell type {cell_type} :",
        is_zi_pred_genelabel_here.mean(),
        "\n",
    )


################
merged= ad.read_h5ad("/beegfs/scratch/ric.boletta/ric.boletta/tsc1/RCC_integrated/python/harmony_P33_P80_integrated_hvg_3000_no_outliers_n10.h5ad")


IC= merged[merged.obs["celltype"].isin(["KO_IC", "KO_transitional_regular"])]

IC =IC[IC.obs.pct_counts_mt <= 30, :]
IC.obs = IC.obs.iloc[:,:30]

IC
#View of AnnData object with n_obs × n_vars = 2514 × 27836

print(IC.obs["celltype"].value_counts())
#ccelltype
#KO_transitional_regular    1476
#KO_IC                      1038
#Name: count, dtype: int64

IC.obs.loc[:,["batch", "celltype"]].groupby("batch").value_counts().reset_index().pivot(index="batch", columns="celltype", values="count").transpose()
#batch                      1    2    3    4
#celltype                                   
#KO_IC                     79  203  235  497
#KO_transitional_regular  183  256  439  578


del IC.layers
IC.layers["counts"]= IC.X.copy()

IC.layers["counts"]= IC.X.copy()

scvi.data.poisson_gene_selection(
    IC,
    n_top_genes=4000,
    batch_key="batch",
    subset=True,
    layer="counts",
)

scvi.model.AUTOZI.setup_anndata(
    IC,
    labels_key="celltype",
    batch_key="batch",
    layer="counts",
)
model = scvi.model.AUTOZI(IC)
model.train(max_epochs=200, plan_kwargs={"lr": 1e-2})

outputs = model.get_alphas_betas()
alpha_posterior = outputs["alpha_posterior"]
beta_posterior = outputs["beta_posterior"]

# Threshold (or Kzinb/Knb+Kzinb in paper)
threshold = 0.5

# q(delta_g < 0.5) probabilities
zi_probs = beta.cdf(0.5, alpha_posterior, beta_posterior)

# ZI genes
is_zi_pred = zi_probs > threshold

print("Fraction of predicted ZI genes :", is_zi_pred.mean())
#Fraction of predicted ZI genes : 0.9395

threshold = 0.75
zi_probs = beta.cdf(0.75, alpha_posterior, beta_posterior)
is_zi_pred = zi_probs > threshold
print("Fraction of predicted ZI genes :", is_zi_pred.mean())
#Fraction of predicted ZI genes : 0.934
threshold = 0.8
#q(delta_g < 0.8) probabilities
zi_probs = beta.cdf(0.8, alpha_posterior, beta_posterior)

is_zi_pred = zi_probs > threshold
print("Fraction of predicted ZI genes :", is_zi_pred.mean())à
Fraction of predicted ZI genes : 0.93375

Analyze gene-cell-type-specific ZI

One may argue that zero-inflation should also be treated on the cell-type (or ‘label’) level, in addition to the gene level. AutoZI can be extended by assuming a random variable 
 for each gene 
 and cell type 
 which denotes the probability that gene 
 is not zero-inflated in cell-type 
. The analysis above can be extended to this new scale.

# Model definition
model_genelabel = scvi.model.AUTOZI(IC, dispersion="gene-label", 
                                    zero_inflation="gene-batch")

# Training
model_genelabel.train(max_epochs=200, plan_kwargs={"lr": 1e-2})

# Retrieve posterior distribution parameters
outputs_genelabel = model_genelabel.get_alphas_betas()
alpha_posterior_genelabel = outputs_genelabel["alpha_posterior"]
beta_posterior_genelabel = outputs_genelabel["beta_posterior"]

threshold = 0.5
# q(delta_g < 0.5) probabilities
zi_probs_genelabel = beta.cdf(0.5, alpha_posterior_genelabel, beta_posterior_genelabel)

# ZI gene-cell-types
is_zi_pred_genelabel = zi_probs_genelabel > threshold

ct = IC.obs.celltype.astype("category")
codes = np.unique(ct.cat.codes)
cats = ct.cat.categories
for ind_cell_type, cell_type in zip(codes, cats, ):
    is_zi_pred_genelabel_here = is_zi_pred_genelabel[:, ind_cell_type]
    print(
        f"Fraction of predicted ZI genes for cell type {cell_type} :",
        is_zi_pred_genelabel_here.mean(),
        "\n",
    )



for ind_cell_type, cell_type in zip(codes, cats):
    mask_sufficient_expression = (
        np.array(IC.X[IC.obs.celltype.values.reshape(-1) == cell_type, :].mean(axis=0))
        > 1.0
    ).reshape(-1)
    print(
        f"Fraction of genes with avg expression > 1 for cell type {cell_type} :",
        mask_sufficient_expression.mean(),
    )
    is_zi_pred_genelabel_here = is_zi_pred_genelabel[mask_sufficient_expression, ind_cell_type]
    print(
        f"Fraction of predicted ZI genes with avg expression > 1 for cell type {cell_type} :",
        is_zi_pred_genelabel_here.mean(),
        "\n",
    )

