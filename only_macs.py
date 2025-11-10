###only macs
#cd ..//..//beegfs/scratch/ric.boletta/donato.alessia/immune_P33_P80/
#srun -p interactive --mem=64GB --x11 --pty bash
#conda activate python39

import os
import anndata as ad
import scanpy as sc
import pandas as pd
import scanpy.external as sce
import schist as scs
import matplotlib as plt
plt.use('tkagg')


cluster_level2= ["Prox_tub",
"eGFP+",
"Activate_mono",
"Tcell",
"Spp1+_MT",
"Spp1_eGFP",
"Tcell_2", 
"Macs_monoder",
"Damage_cell",
"Macs",
"Neutrophils", 
"Dc_monoder",
"Macs_MRC1"]

macs= immune[immune.obs["celltype_annotation"].isin(["Activate_mono", "Macs_monoder", "Macs", "Macs_MRC1"])]
macs.X= macs.layers["counts"]
macs.obs = macs.obs.iloc[:,:14]
del macs.layers
sc.pl.scatter(macs, x="total_counts", y="pct_counts_mt")
sc.pl.scatter(macs, x="total_counts", y="n_genes_by_counts")

macs.raw = macs
macs.layers["counts"] = macs.X
sc.pl.scatter(macs, x="total_counts", y="pct_counts_mt")

sc.pp.normalize_total(macs)
macs.layers["cpm"] = macs.X
sc.pp.log1p(macs)
macs.layers["log"] = macs.X

sc.pp.highly_variable_genes(macs, n_top_genes =3000)

sc.pp.pca(macs, n_comps=50, use_highly_variable=True)
sc.pp.scale(macs)
macs.layers["scaled"] = macs.X
sc.tl.pca(macs, svd_solver='arpack')
sc.pl.pca_variance_ratio(macs, log=True, n_pcs=50) 

sce.pp.harmony_integrate(macs, key="batch")

macs.obsm['X_pca'] = macs.obsm['X_pca_harmony']

sc.pp.neighbors(macs, n_neighbors=10, n_pcs=25)
sc.tl.umap(macs)
sc.pl.umap(macs, color="batch", save="macs_batch.png")

sc.tl.leiden(macs, key_added="leiden_res0_25", resolution=0.25)
sc.tl.leiden(macs, key_added="leiden_res0_5", resolution=0.5)
sc.tl.leiden(macs, key_added="leiden_res1", resolution=1.0)
sc.tl.leiden(macs, key_added="leiden_res1_25", resolution=1.25)
sc.tl.leiden(macs, key_added="leiden_res0_4", resolution=0.4)

sc.pl.umap(macs, color=['leiden_res0_25', 'leiden_res0_4','leiden_res0_5', 'leiden_res1', 'leiden_res1_25'], 
           ncols=2, legend_loc='on data')

sc.pl.umap(macs, color=["leiden_res0_25", "Il1b", "Il1a"])
scs.inference.fit_model(macs, n_init=100)
sc.pl.umap(macs, color=['nsbm_level_0', 'nsbm_level_1', 'nsbm_level_2', 'nsbm_level_3', 'nsbm_level_4'], 
           ncols=2, legend_loc='on data')

sc.tl.rank_genes_groups(macs, layer = "counts",use_raw=False,  method="wilcoxon", groupby="leiden_res0_25")
pd.DataFrame(macs.uns["rank_genes_groups"]["names"]).head(15)
         0       1        2        3       4
0     C1qb    Atf3   S100a6    Fxyd2  Hspa1b
1     C1qa    Cd83      Vim     Gpx3  Hspa1a
2     Apoe    Il1b    Crip1    Cox6c    Rpl8
3    Itm2b  Nfkbiz   Tmsb10     Spp1    Lyz2
4     Lyz2    Cd14   S100a4      Dbi     Pf4
5     Ctsb   Dusp1  S100a11    Cox5a    Ptma
6     Lgmn  Nfe2l2    Ahnak   Ndufa4  Ifitm2
7   Tmsb4x    Pim1     Lsp1    Cox8a     Jun
8    Rps29  Nfkbia     Ccr2    Cox7c    Sem1
9    Txnip   Zfp36  S100a10  Chchd10    Lifr
10  Hspa1a    Btg2     Rpsa    Cox5b    Oaz1
11    C1qc     Fos    Napsa   Atp5g1  Ifitm3
12  Fcer1g    Rgs1    Anxa2   Cox7a2  Arpp19
13    Aif1    Junb    Rplp0   Uqcr10    H1f5
14   Ms4a7    Ccl4     Gm2a     Gpx4   Smdt1

sc.tl.rank_genes_groups(macs, layer = "counts",use_raw=False,  method="wilcoxon", groupby="leiden_res0_5")
pd.DataFrame(macs.uns["rank_genes_groups"]["names"]).head(15)
         0        1       2        3       4        5       6         7        8         9
0     C1qb   Nfe2l2    Junb   S100a6     Kap     Ctsd   Fxyd2     Ifit2  S100a11    Hspa1b
1     C1qa     Atf3   Socs3    Crip1   H2-Aa     Ftl1    Spp1     Isg15   S100a8    Hspa1a
2    Itm2b     Il1b     Fos      Vim    Cd74     Ctsb    Gpx3     Ifit3    Prdx6      Rpl8
3     Apoe     Cd83     Jun   S100a4  H2-Ab1     Cd63   Cox6c    Ifi204   Vps37b       Jun
4   Tmsb4x     Rgs1    Jund     Ccr2  H2-Eb1    Trem2     Dbi     Slfn5   S100a6      Lyz2
5     Ctsb   Neurl3   Dusp1   Tmsb10    Ccl4     Fth1  Ndufa4    Ifitm3   S100a9  Hsp90aa1
6     Lyz2   Nfkbiz    Ier3    Ahnak   Dusp1     Lyz2   Cox7c      Irf7    Rps27    Ifitm3
7     Aif1  Bcl2a1b  Nfkbia     Lsp1    Ccl3     Cstb   Cox8a    Rnf213     Eif1       Pf4
8    Ms4a7     Btg2   Zfp36     Rpsa   S100g  Atp6v0c   Cox5a  Lgals3bp     Btg1    Ifitm2
9   Tyrobp     Cd14     Mt1  S100a10  Lilra5    Fcrls  Cox7a2     Mndal     Ccl5      Rps3
10    C1qc   Tent5c    Ier2  S100a11    Klk1    Pmp22  mt-Nd5    Parp14     Txn1      Lifr
11    Cd52     Fosb   H3f3b    Napsa    Klf2      Grn    Gpx4      Bst2     Cd3g      Jund
12    Cyba     Btg1  Hspa1a    Rplp0   H2-K1     Igf1  Atp5g1     Oasl2   Tmsb10     H3f3b
13  Fcer1g     Skil    Egr1    Anxa2   Fcgr4    Lamp1   Cox5b      Ly6e  Slc7a11      Klf6
14  Man2b1  Tnfaip3    Atf3     Tpt1     Ubc     Ctsz  Uqcr10    Ms4a4c    Dusp5    Arpp19



#### uso raw counts

macs.X = macs.layer["counts"]

sc.pl.umap(macs, color=["leiden_res0_25","leiden_res0_5", "Il1b", "Il1a"], ncols=2)

sc.tl.rank_genes_groups(macs, layer = "counts",use_raw=False,  method="wilcoxon", groupby="leiden_res0_5")
pd.DataFrame(macs.uns["rank_genes_groups"]["names"]).head(15)

 0       1       2       3       4        5        6        7       8        9        10       11
0     Apoe    Rgs1   H2-Aa   Socs3    Spp1   S100a6   mt-Nd2     Gpx3    Spp1      Kap     Slfn5     Lsp1
1     C1qb  Nfe2l2    Cd74  Nfkbia    Ctsd    Crip1   H2-Eb1    Fxyd2    Gpx3    S100g     Ifit2    Psma7
2     C1qa    Atf3  H2-Ab1     Fau    Ftl1   S100a4   mt-Nd4  S100a11   F13a1     Klk1      Irf7   Tbc1d4
3     Lyz2    Fosb  H2-Eb1   Rps24    Cd63      Vim  mt-Cytb    Aldob   Cox6c  mt-Cytb      Ctss    Itga4
4    Itm2b    Il1b     Kap    Junb     Mif     Ccr2    H2-Aa     Mpc2   Mkrn1   mt-Nd2     Isg15    Crip1
5     Ctsb  Neurl3   Dusp1    Rps9   Uqcrq    Ahnak   H2-Ab1      Dbi    Gpx4   Cited2     Ifit3  Ccdc88a
6    Rps29    Btg2    Ccl4    Rpl8   Cox6c   Tmsb10     Cd74     Gpx4   Fxyd2   mt-Nd4  Lgals3bp  Tbc1d15
7   Tmsb4x  Tent5c    Ccl3     Fos   Trem2     Lsp1   mt-Co3    Prdx6    Ier3    H2-Aa    Parp14    H2-Q6
8     Lgmn    Cd83  Lilra5   Rpl13  Cox7a2     Rpsa   mt-Nd1    Cox6c   Rps27   mt-Nd5      Rgs1      Vim
9    Trem2  Pmepa1     Ubc   Zfp36   Aldoa  S100a10   mt-Co2   Cox7a2    Sdc4   mt-Nd1      Ly6e    Cox5a
10    Cyba    Rgs2   Mmp13   Rpl32    Gpx4  S100a11  mt-Atp6   Vps37b  Hspa1a     Cd74    Ifi204  Gadd45b
11  Tyrobp    Junb   Fcgr4    Rps7    Fth1    Napsa   mt-Co1   Tmsb10  Atp5a1   H2-Eb1      Rgs2    Myo1g
12   Rps23    Btg1   S100g   Rps18   Prdx1    Rplp0   mt-Nd3   Ndufa4   Cox7c   mt-Nd3      Psap    Socs2
13   Rpl11    Ctss   H2-K1   Rps13  Uqcr10     Tpt1   mt-Nd5   S100a6     Mt1     Ccl4     Mndal   Tmsb10
14    Cd63    Skil  mt-Nd2    Jund  Atp5g3    Anxa2  H2-DMb1   S100a8    Klf9   H2-Ab1      Apoe    Fabp5



sc.pl.rank_genes_groups_dotplot(macs, groupby="leiden_res0_5", 
                                standard_scale="var", n_genes=5, save="macs_top5.png")

