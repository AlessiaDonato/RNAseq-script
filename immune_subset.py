##subset immune cell from integrated data
##
#cd ..//..//beegfs/scratch/ric.boletta/donato.alessia/immune_P33_P80/
#srun -p interactive --mem=64GB --x11 --pty bash
#conda activate python39


import os
import anndata as ad
import scanpy as sc
import pandas as pd
import scanpy.external as sce
import schist as scs
import matplotlib
matplotlib.use('tkagg')

merged= ad.read_h5ad("/beegfs/scratch/ric.boletta/ric.boletta/tsc1/RCC_integrated/python/harmony_P33_P80_integrated_hvg_3000_no_outliers_n10.h5ad")

immune= merged[merged.obs["celltype"].isin(["Macrophages", "Neutrophiles", "T_CD4", "T_CD8", "B", "T", "NK"])]
immune.obs = immune.obs.iloc[:,:14]

sc.pl.scatter(immune, x="total_counts", y="pct_counts_mt")
sc.pl.scatter(immune, x="total_counts", y="n_genes_by_counts")
immune.raw = immune
del immune.layers
immune = immune[immune.obs.pct_counts_mt < 30, :]
#check if X is row counts 
immune.layers["counts"] = immune.X.copy
sc.pl.scatter(immune, x="total_counts", y="pct_counts_mt")

sc.pp.normalize_total(immune)
immune.layers["cpm"] = immune.X.copy()
sc.pp.log1p(immune)
immune.layers["log"] = immune.X.copy()

sc.pp.highly_variable_genes(immune)


sc.pp.pca(immune, n_comps=50, use_highly_variable=True)

#immune = immune[:,immune.var.highly_variable]
sc.pp.scale(immune)
immune.layers["scaled"] = immune.X

sc.tl.pca(immune, svd_solver='arpack')
sc.pl.pca_variance_ratio(immune, log=True, n_pcs=50) 

sce.pp.harmony_integrate(immune, key="batch")

immune.obsm['X_pca'] = immune.obsm['X_pca_harmony']
sc.pp.neighbors(immune, n_neighbors=10, n_pcs=40)
sc.tl.umap(immune)
sc.pl.umap(immune, color="batch")
sc.tl.leiden(immune, resolution=0.5)

scs.inference.fit_model(immune, n_init=100)
immune.write_h5ad("immune_harmony_umap.h5ad")

sc.pl.umap(immune, color=['nsbm_level_0', 'nsbm_level_1', 'nsbm_level_2', 'nsbm_level_3', 'nsbm_level_4'], 
           ncols=2, legend_loc='on data')

sc.tl.paga(immune, groups='nsbm_level_2')
sc.pl.paga(immune, color = 'nsbm_level_2')
###umap 3d
##salvo la vecchia umap 
immune.obsm["umap_initial"] = immune.obsm["X_umap"]

### calcolo pca con 3 componenti
sc.tl.umap(immune, n_components=3)
immune.obsm["X_umap"]
array([[-1.4212736 , 12.865025  , -5.719051  ],
       [ 4.7941027 ,  4.608547  ,  3.6009204 ],
       [-0.0153746 ,  4.218623  ,  1.3205414 ],
       ...,
       [ 0.5809965 ,  3.765774  , -1.5452048 ],
       [ 0.9622101 ,  1.5386674 , -0.46414635],
       [-0.37935337,  8.246139  , 10.474813  ]], dtype=float32)

sc.pl.embedding(immune, basis="X_umap", projection="3d", color="nsbm_level_2", s=25 )

####Gene marker
##use counts

sc.tl.rank_genes_groups(immune, layer = immune.layers["counts"], method="wilcoxon", groupby="nsbm_level_2")

pd.DataFrame(immune.uns["rank_genes_groups"]["names"]).head(15)
          0        1       2       3        4        5        6        7        8       9       10       11       12
0     Aldob    Wfdc2    Gm2a  Tmsb10     Spp1   S100a9     Cd3g    Trem2    Sparc    Cd74   S100a8     Lyz2     Apoe
1    Spink1     Ldhb  Wfdc17  Ms4a4b   mt-Nd5   S100a8     Cd3d     Apoe    Fabp4    C1qa   S100a9      Vim     Mrc1
2   Slc34a1     Spp1  Tmsb4x   Trbc2    Cox6c    Il1r2   Tmsb10      Pf4    Plvap    C1qb    Il1r2   Lgals3     C1qc
3      Gpx3   Atp1b1  Laptm5    Cd3g    Fxyd2    Fxyd2   S100a4     Aif1    Gng11    C1qc      Hdc     Emp3     C1qa
4   Mir6236  Chchd10  Lgals3    Nkg7    Wfdc2     Spp1     Cd3e     Lyz2   Hba-a1    Ctss     Mxd1   Tyrobp     Lyz2
5       Kap    Cox6c    Ctss    Ccl5   mt-Nd1  S100a11  Ptprcap     C1qb     Flt1    Apoe       Hp   Fcer1g     C1qb
6     Ttc36    mEGFP  H2-DMa    Cd3d   Chchd2     Ldhb  S100a10    Blvrb   Col4a1  H2-Eb1  S100a11     Cyba      Trf
7     Acsm2   S100a1    Cd52  Gimap4   mt-Nd4     Lcn2    Cxcr6     C1qa    Egfl7   H2-Aa     Srgn    Crip1    Vcam1
8    mt-Nd2   Atp5g3    Cst3  Gimap3   mt-Co2    Wfdc2    Hmgb1   Cxcl14   Hbb-bs  H2-Ab1   Clec4d     Cybb    Csf1r
9     Fxyd2   Atp1a1    Rps9  Gimap6    Cox7c    mEGFP     Trdc    Gins2   Hbb-bt   Ms4a7   S100a6   S100a4   Ms4a6c
10  mt-Cytb   Ndufa4    Ctsz  Vps37b  mt-Cytb   S100a6      Vim    Cenps   Hba-a2    Cybb   Retnlg    Mpeg1      Pf4
11    Cltrn   Ndufb9   Crip1    Ets1    Uqcrb   Retnlg      Ltb    Fcrls  Gpihbp1    Rgs1     Slpi   Wfdc17  Selenop
12    Chpt1    Cd24a  Tyrobp    Hcst  mt-Atp6  Tnfaip2     Lsp1    Rad51   Igfbp3    Ctsh    Cxcr2  Alox5ap   Pla2g7
13     Acy3    Cox8a    Actb   H2-Q7   mt-Nd2  Chchd10     Ets1  Psmc3ip   Adgrl4    Mafb    Csf3r      Msn     Ctsb
14     Lrp2    Cox5a     Vim    Trac    Cox5b      Hdc    Crip1    Brip1    Lars2   Csf1r     G0s2     Zeb2    Fcgr3


sc.pl.dotplot(immune, "Ptprc", groupby="nsbm_level_2", save="dotplot_CD45.png")
sc.pl.umap(immune,color="nsbm_level_2",palette= "tab20b", legend_loc="on data", save="Umap_immune.png")
print(immune.obs["nsbm_level_2"].value_counts())
nsbm_level_2
9     2730
10    1392
3     1344
12    1097
1      977
0      681
11     565
2      491
7      158
4      150
6       84
8       60
5       30

## creo lista marker per identificazione cluster
##MNPs, Tcell, Bcell, neutrophils, endo, S1/S2proximal tubule, S3 proximal tubule, interstitial cells
marker= { "MNPs" : ["Lyz2", "C1qa"] , "Tcell":["Ccl5", "Cd3g"], "Bcell": ["Igkc", "Cd79a"], "neutrophils":["S100a9", "S100a8"], "other":["Emcn", "Meis2", "Miox", "Aldob", "Kap", "Napsa", "Acta2"]}

sc.pl.dotplot(immune, marker, groupby="nsbm_level_2", save="dotplot_marker.png")

endo= ["Emcn", "Meis2"]
sc.pl.dotplot(immune, endo, groupby="nsbm_level_2")

marker_2= {"Macs" : ["C1qa", "C1qb", "C1qc","Cd14", "Cd72"], "KMRs": ["Cx3cr1", "Cd81","Adgre1", "Itgam", "C1qa"], "Ly6c_hi":["Plac8", "Chil3", "S100a4", "S100a6", "Hp"], "Ly6c_lo": ["Ear2", "Pglyrp1", "Cebpb1", "Itgal1", "Eno3"], "Mrc1+ KRMs":["Mrc1","Ccl8","Fcrls","Pf4", "C4b", "Ccl7"],
           "Spp1+ KRMs": ["Gm42418", "Gpx3","mt-Co3", "mt-Nd2", "Ndufb9", "Spp1"], "mono/DC_mixture": ["Actg1","Crip1", "Hspa8", "Rplp0", "Vim"], "cDC2": ["Cd209a", "Clec10","Klrd1","H2-Oa","H2-DMb2"], "pDCs":["Ly6d", "Siglech", "Cox6a2", "Ccr9", "Iglc3"], "CD8+ T": ["Cd8a", "Gzmk"], "Tregs": ["Ikzf2", "Foxp3"], 
"effector CD4+": ["Cd4", "Cxcr3", "Cd40lg"], "central memory T":["Sell", "Ccr7", "Lef1", "Cd69"], "CD4+ Th17": ["Rorc", "Tmem176a", "Lgals3"], "NKT1": ["Klrk1", "Ly6c2", "Anxa1"], 
"Gzma+ NK":["Gzma", "Tyrobp"], "Gzma+ CD8+":["Gzma", "Cx3cr1", "Zeb2"], "Gzma_lo NK":["Cd7", "Xcl1"], "CD4+ T cells_IFNγ": ["Ifit1", "Ifit3"], "ILCs":["Il1rl1", "Arg1"], 
"B1 B cells": ["Ccr7", "Cd69"], "T3/follicular helper B cells":["Cr2", "Ccr6"], "T1 B cells/memory B cells": ["Spib", "Iglc1"], "plasma cells": ["Zbtb20", "Apoe"], 
"memory B cells": ["Ifit3", "Isg15"]}

#KeyError: "Could not find keys '['Cebpb1', 'Clec10', 'Gm42418', 'Itgal1']' in columns of `adata.obs` or in adata.raw.var_names."

###controllo che tutti i marker siano nei geni del dataset
filtered_markers = {key: [gene for gene in values if gene in immune.var_names] for key,values in marker_2.items() }
sc.pl.dotplot(immune, filtered_markers, groupby="nsbm_level_2")
sc.pl.dotplot(immune, filtered_markers, groupby=["nsbm_level_2", "condition", "timepoints"])
sc.pl.dotplot(immune, "Lox", groupby="nsbm_level_2")
sc.pl.dotplot(immune, filtered_markers, groupby="nsbm_level_2", layer="scaled")


fig, axes = plt.subplots(1,2)
sc.pl.umap(immune,color="nsbm_level_2", s=25, ax=axes[0],show=False, legend_loc="on data")
sc.pl.rank_genes_groups(immune,ax=axes[1])


sc.pl.umap(immune, color=['nsbm_level_0', 'nsbm_level_1', 'nsbm_level_2', 'nsbm_level_3', 'nsbm_level_4'], 
           ncols=2, legend_loc='on data')

sc.pl.umap(immune, color='nsbm_level_0')

macs= {"Macs" : ["C1qa", "C1qb", "C1qc","Cd14", "Cd72"], "KMRs": ["Cx3cr1", "Cd81", "Adgre1", "Itgam", "C1qa"], "Ly6c_hi":["Plac8", "Chil3", "S100a4", "S100a6", "Hp"],
"Ly6c_lo": ["Ear2", "Pglyrp1", "Cebpb1", "Itgal1", "Eno3"],"Mrc1+ KRMs":["Mrc1", "Ccl8","Fcrls","Pf4", "C4b", "Ccl7"],"Spp1+ KRMs": ["Gm42418", "Gpx3","mt-Co3", "mt-Nd2", "Ndufb9", "Spp1"]}
filtered_macs = {key: [gene for gene in values if gene in immune.var_names] for key,values in macs.items() }


# In this example we want to show UMAPs of different cell type markers,
# with markers of a single cell type in one row
# and with a different number of markers per cell type (row)

import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.colors as mcolors

# Inital setting for plot size
from matplotlib import rcParams

FIGSIZE = (3, 3)
rcParams["figure.figsize"] = FIGSIZE

# Make Axes
# Number of needed rows and columns (based on the row with the most columns)
nrow = len(filtered_macs)
ncol = max([len(vs) for vs in filtered_macs.values()])
fig, axs = plt.subplots(nrow, ncol, figsize=(2 * ncol, 2 * nrow))
# Plot expression for every marker on the corresponding Axes object
for row_idx, (cell_type, markers) in enumerate(filtered_macs.items()):
    col_idx = 0
    for marker in markers:
        ax = axs[row_idx, col_idx]
        sc.pl.umap(immune, color=marker, ax=ax, show=False, frameon=False, s=20)
        # Add cell type as row label - here we simply add it as ylabel of
        # the first Axes object in the row
        if col_idx == 0:
            # We disabled axis drawing in UMAP to have plots without background and border
            # so we need to re-enable axis to plot the ylabel
            ax.axis("on")
            ax.tick_params(
                top="off",
                bottom="off",
                left="off",
                right="off",
                labelleft="on",
                labelbottom="off",
            )
            ax.set_ylabel(cell_type + "\n", rotation=90, fontsize=14)
            ax.set(frame_on=False)
        col_idx += 1
    # Remove unused column Axes in the current row
    while col_idx < ncol:
        axs[row_idx, col_idx].remove()
        col_idx += 1
# Alignment within the Figure
fig.tight_layout()
sc.pl.umap(immune)

###########DENDRITIC cell##########

DCs= {"cDC2": ["Cd209a", "Clec10","Klrd1","H2-Oa","H2-DMb2"],"pDCs":["Ly6d", "Siglech", "Cox6a2", "Ccr9", "Iglc3"], "mono/DC_mixture": ["Actg1","Crip1", "Hspa8", "Rplp0", "Vim"]}
filtered_DCs = {key: [gene for gene in values if gene in immune.var_names] for key,values in DCs.items() }

# In this example we want to show UMAPs of different cell type markers,
# with markers of a single cell type in one row
# and with a different number of markers per cell type (row)

import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.colors as mcolors

# Inital setting for plot size
from matplotlib import rcParams

FIGSIZE = (4, 4)
rcParams["figure.figsize"] = FIGSIZE

# Make Axes
# Number of needed rows and columns (based on the row with the most columns)
nrow = len(filtered_DCs)
ncol = max([len(vs) for vs in filtered_DCs.values()])
fig, axs = plt.subplots(nrow, ncol, figsize=(2 * ncol, 2 * nrow))
# Plot expression for every marker on the corresponding Axes object
for row_idx, (cell_type, markers) in enumerate(filtered_DCs.items()):
    col_idx = 0
    for marker in markers:
        ax = axs[row_idx, col_idx]
        sc.pl.umap(immune, color=marker, ax=ax, show=False, frameon=False, s=20)
        # Add cell type as row label - here we simply add it as ylabel of
        # the first Axes object in the row
        if col_idx == 0:
            # We disabled axis drawing in UMAP to have plots without background and border
            # so we need to re-enable axis to plot the ylabel
            ax.axis("on")
            ax.tick_params(
                top="off",
                bottom="off",
                left="off",
                right="off",
                labelleft="on",
                labelbottom="off",
            )
            ax.set_ylabel(cell_type + "\n", rotation=90, fontsize=14)
            ax.set(frame_on=False)
        col_idx += 1
    # Remove unused column Axes in the current row
    while col_idx < ncol:
        axs[row_idx, col_idx].remove()
        col_idx += 1
# Alignment within the Figure
fig.tight_layout()
sc.pl.umap(immune)


#############################################################
marker_T= {"CD8+ T": ["Cd8a", "Gzmk"], "Tregs": ["Ikzf2", "Foxp3"], "effector CD4+": ["Cxcr3", "Cd40lg"], "central memory T":["Sell", "Ccr7", "Lef1", "Cd69"], 
"CD4+ Th17": ["Rorc", "Tmem176a", "Lgals3"],"CD4+ T cells_IFNγ": ["Ifit1", "Ifit3"], "Gzma+ CD8+":["Gzma", "Cx3cr1", "Zeb2"]}

filtered_T = {key: [gene for gene in values if gene in immune.var_names] for key,values in marker_T.items() }

 # In this example we want to show UMAPs of different cell type markers,
# with markers of a single cell type in one row
# and with a different number of markers per cell type (row)

import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.colors as mcolors

# Inital setting for plot size
from matplotlib import rcParams

FIGSIZE = (10, 10)
rcParams["figure.figsize"] = FIGSIZE

# Make Axes
# Number of needed rows and columns (based on the row with the most columns)
nrow = len(filtered_T)
ncol = max([len(vs) for vs in filtered_T.values()])
fig, axs = plt.subplots(nrow, ncol, figsize=(2 * ncol, 2 * nrow))
# Plot expression for every marker on the corresponding Axes object
for row_idx, (cell_type, markers) in enumerate(filtered_T.items()):
    col_idx = 0
    for marker in markers:
        ax = axs[row_idx, col_idx]
        sc.pl.umap(immune, color=marker, ax=ax, show=False, frameon=False, s=20)
        # Add cell type as row label - here we simply add it as ylabel of
        # the first Axes object in the row
        if col_idx == 0:
            # We disabled axis drawing in UMAP to have plots without background and border
            # so we need to re-enable axis to plot the ylabel
            ax.axis("on")
            ax.tick_params(
                top="off",
                bottom="off",
                left="off",
                right="off",
                labelleft="on",
                labelbottom="off",
            )
            ax.set_ylabel(cell_type + "\n", rotation=90, fontsize=14)
            ax.set(frame_on=False)
        col_idx += 1
    # Remove unused column Axes in the current row
    while col_idx < ncol:
        axs[row_idx, col_idx].remove()
        col_idx += 1
# Alignment within the Figure
fig.tight_layout()
sc.pl.umap(immune)

marker_NK={"NKT1": ["Klrk1", "Ly6c2", "Anxa1"], "Gzma+ NK":["Gzma", "Tyrobp"], "Gzma_low NK":["Cd7", "Xcl1"], "ILCs":["Il1rl1", "Arg1"]}
filtered_NK = {key: [gene for gene in values if gene in immune.var_names] for key,values in marker_NK.items() }

tmp=pd.crosstab(immune.obs['nsbm_level_2'],immune.obs['condition'], normalize='columns').T.plot(kind='bar', stacked=True)
tmp.legend(title='nsbm_level_2', bbox_to_anchor=(1.26, 1.02),loc='upper right')

cell_proportion_df = pd.crosstab(immune.obs['nsbm_level_2'],immune.obs['condition'], normalize='columns').T.plot(kind='bar', stacked=True)

##########
immune.obs.loc[:,["nsbm_level_2","pct_counts_mt"]]
immune.obs.loc[:,["nsbm_level_2","pct_counts_mt"]].groupby("nsbm_level_2").mean()
####tolgo geni mitocondriale dai più variabili
immune.var.loc[[x.startswith("mt-") for x in immune.var_names] ,"highly_variable"] = False

sc.pp.pca(immune, n_comps=50, use_highly_variable=True)


sc.pp.scale(immune)
immune.layers["scaled"] = immune.X

sc.tl.pca(immune, svd_solver='arpack')
sc.pl.pca_variance_ratio(immune, log=True, n_pcs=50) 

sce.pp.harmony_integrate(immune, key="batch")

immune.obsm['X_pca'] = immune.obsm['X_pca_harmony']
sc.pp.neighbors(immune, n_neighbors=10, n_pcs=35)
sc.tl.umap(immune)
sc.pl.umap(immune, color="batch")
sc.tl.leiden(immune, resolution=1, key_added = "leiden_1")
sc.tl.leiden(immune, resolution=0.8, key_added = "leiden_0.8")
sc.tl.leiden(immune, resolution=0.4, key_added = "leiden_0.4")
sc.tl.leiden(immune, resolution=1.4, key_added = "leiden_1.4")
sc.pl.umap(immune, color=["leiden_1.4", "leiden_1", "leiden_0.8", "leiden_0.4"])


scs.inference.fit_model(immune, n_init=100)
#immune.write_h5ad("immune_harmony_umap.h5ad")

sc.pl.umap(immune_no_mito_genes, color=['nsbm_level_0', 'nsbm_level_1', 'nsbm_level_2', 'nsbm_level_3', 'nsbm_level_4'], 
           ncols=2, legend_loc='on data')
sc.pl.umap(immune_no_mito_genes, color=['nsbm_level_0', 'nsbm_level_1', 'nsbm_level_2', 'nsbm_level_3', 'nsbm_level_4'], 
           ncols=2, legend_loc='on data')


sc.tl.rank_genes_groups(immune, use_raw=True, method="wilcoxon", groupby="leiden_0.4")

pd.DataFrame(immune.uns["rank_genes_groups"]["names"]).head(10)

sc.pl.dotplot(immune, filtered_markers, groupby="leiden_0.4", layer="scaled")
#per vedere qunante cellule in ogni cluster vengono dai vari batch
immune.obs.loc[:, ["leiden_0.4","batch_id"]].groupby("batch_id").value_counts().reset_index().pivot(columns="batch_id",index="leiden_0.4")
#percentuale mt in ogni cluster  per batch
immune.obs.loc[:, ["leiden_0.4","batch_id","pct_counts_mt"]].groupby(["leiden_0.4","batch_id"]).mean().reset_index().pivot(columns="batch_id", index="leiden_0.4")


###umap 3d
##salvo la vecchia umap 
immune.obsm["umap_initial"] = immune.obsm["X_umap"]
### calcolo pca con 3 componenti
sc.tl.umap(immune, n_components=3)
immune.obsm["X_umap"]
array([[-1.4212736 , 12.865025  , -5.719051  ],
       [ 4.7941027 ,  4.608547  ,  3.6009204 ],
       [-0.0153746 ,  4.218623  ,  1.3205414 ],
       ...,
       [ 0.5809965 ,  3.765774  , -1.5452048 ],
       [ 0.9622101 ,  1.5386674 , -0.46414635],
       [-0.37935337,  8.246139  , 10.474813  ]], dtype=float32)

sc.pl.embedding(immune, basis="X_umap", projection="3d", color="nsbm_level_2", s=25 )


