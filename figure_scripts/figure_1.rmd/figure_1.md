---
title: "Figure 1"
output: rmarkdown::github_document
---


```r
library(reticulate)
use_python("/home/nealpsmith/.conda/envs/sc_analysis/bin/python")
```

First, we need to aggregate the matrices and assess the data quality.  Here we use Pegasus to create a count matrix and to
calcluate some quality-control statistics (% mitochondrial UMIs, number of genes per cell).  From there, we can plot
these statistics to determine the proper cutoffs to define the quality cells.  Based on the distributions, we chose to include
cells with < 30% mitochondrial reads and > 500 genes. test.


```python
import pegasus as pg
import scanpy as sc
import pandas as pd
import matplotlib.colors as clr
import matplotlib.pyplot as plt
import matplotlib as mpl

# Set a colormap
colormap = clr.LinearSegmentedColormap.from_list('gene_cmap', ["#d3d3d3" ,'#482cc7'], N=200)

# Aggregate the matrices
pg.aggregate_matrices(csv_file = "/home/nealpsmith/projects/medoff/cellranger/aggregate_matrix.csv",
                      what_to_return = "/home/nealpsmith/projects/medoff/data/all_data.h5sc")
```

```
## 2021-10-01 10:06:17,898 - pegasus - INFO - Time spent on 'read_input' = 0.22s.
## Processed /home/nealpsmith/projects/medoff/data/h5_files/500008_ANA_Ag_filtered_feature_bc_matrix.h5.
## 2021-10-01 10:06:18,068 - pegasus - INFO - Time spent on 'read_input' = 0.17s.
## Processed /home/nealpsmith/projects/medoff/data/h5_files/500008_ANA_Dil_filtered_feature_bc_matrix.h5.
## 2021-10-01 10:06:18,253 - pegasus - INFO - Time spent on 'read_input' = 0.18s.
## Processed /home/nealpsmith/projects/medoff/data/h5_files/500008_ANA_Pre_filtered_feature_bc_matrix.h5.
## 2021-10-01 10:06:18,552 - pegasus - INFO - Time spent on 'read_input' = 0.30s.
## Processed /home/nealpsmith/projects/medoff/data/h5_files/500012_ANA_Ag_filtered_feature_bc_matrix.h5.
## 2021-10-01 10:06:18,892 - pegasus - INFO - Time spent on 'read_input' = 0.34s.
## Processed /home/nealpsmith/projects/medoff/data/h5_files/500012_ANA_Dil_filtered_feature_bc_matrix.h5.
## 2021-10-01 10:06:19,268 - pegasus - INFO - Time spent on 'read_input' = 0.37s.
## Processed /home/nealpsmith/projects/medoff/data/h5_files/500012_ANA_Pre_filtered_feature_bc_matrix.h5.
## 2021-10-01 10:06:19,719 - pegasus - INFO - Time spent on 'read_input' = 0.45s.
## Processed /home/nealpsmith/projects/medoff/data/h5_files/500021_AA_Ag_filtered_feature_bc_matrix.h5.
## 2021-10-01 10:06:20,187 - pegasus - INFO - Time spent on 'read_input' = 0.47s.
## Processed /home/nealpsmith/projects/medoff/data/h5_files/500021_AA_Dil_filtered_feature_bc_matrix.h5.
## 2021-10-01 10:06:20,434 - pegasus - INFO - Time spent on 'read_input' = 0.24s.
## Processed /home/nealpsmith/projects/medoff/data/h5_files/500021_AA_Pre_filtered_feature_bc_matrix.h5.
## 2021-10-01 10:06:20,683 - pegasus - INFO - Time spent on 'read_input' = 0.25s.
## Processed /home/nealpsmith/projects/medoff/data/h5_files/500024_ANA_Ag_filtered_feature_bc_matrix.h5.
## 2021-10-01 10:06:20,967 - pegasus - INFO - Time spent on 'read_input' = 0.28s.
## Processed /home/nealpsmith/projects/medoff/data/h5_files/500024_ANA_Dil_filtered_feature_bc_matrix.h5.
## 2021-10-01 10:06:21,339 - pegasus - INFO - Time spent on 'read_input' = 0.37s.
## Processed /home/nealpsmith/projects/medoff/data/h5_files/500024_ANA_Pre_filtered_feature_bc_matrix.h5.
## 2021-10-01 10:06:21,483 - pegasus - INFO - Time spent on 'read_input' = 0.14s.
## Processed /home/nealpsmith/projects/medoff/data/h5_files/500015_ANA_Ag_filtered_feature_bc_matrix.h5.
## 2021-10-01 10:06:21,723 - pegasus - INFO - Time spent on 'read_input' = 0.24s.
## Processed /home/nealpsmith/projects/medoff/data/h5_files/500015_ANA_Pre_filtered_feature_bc_matrix.h5.
## 2021-10-01 10:06:22,203 - pegasus - INFO - Time spent on 'read_input' = 0.48s.
## Processed /home/nealpsmith/projects/medoff/data/h5_files/500030_AA_Ag_filtered_feature_bc_matrix.h5.
## 2021-10-01 10:06:22,446 - pegasus - INFO - Time spent on 'read_input' = 0.24s.
## Processed /home/nealpsmith/projects/medoff/data/h5_files/500030_AA_Dil_filtered_feature_bc_matrix.h5.
## 2021-10-01 10:06:22,622 - pegasus - INFO - Time spent on 'read_input' = 0.17s.
## Processed /home/nealpsmith/projects/medoff/data/h5_files/500030_AA_Pre_filtered_feature_bc_matrix.h5.
## 2021-10-01 10:06:22,950 - pegasus - INFO - Time spent on 'read_input' = 0.33s.
## Processed /home/nealpsmith/projects/medoff/data/h5_files/500032_AA_Ag_filtered_feature_bc_matrix.h5.
## 2021-10-01 10:06:23,358 - pegasus - INFO - Time spent on 'read_input' = 0.41s.
## Processed /home/nealpsmith/projects/medoff/data/h5_files/500032_AA_Pre_filtered_feature_bc_matrix.h5.
## 2021-10-01 10:06:23,770 - pegasus - INFO - Time spent on 'read_input' = 0.41s.
## Processed /home/nealpsmith/projects/medoff/data/h5_files/500035_AA_Pre_filtered_feature_bc_matrix.h5.
## 2021-10-01 10:06:24,010 - pegasus - INFO - Time spent on 'read_input' = 0.24s.
## Processed /home/nealpsmith/projects/medoff/data/h5_files/500035_AA_Ag_filtered_feature_bc_matrix.h5.
## 2021-10-01 10:06:24,946 - pegasus - INFO - Time spent on 'aggregate' = 0.93s.
## 2021-10-01 10:06:32,407 - pegasus - INFO - Time spent on 'write_output' = 7.46s.
## Aggregated 21 files.
```

```python
all_data = pg.read_input("/home/nealpsmith/projects/medoff/data/all_data.h5sc")
```

```
## 2021-10-01 10:06:35,270 - pegasus - INFO - Time spent on 'read_input' = 2.84s.
```

```python
pg.qc_metrics(all_data, percent_mito = 30)

# Plot the percent mito/n genes
fig, ax = plt.subplots(1)
x = all_data.obs["n_genes"]
y = all_data.obs["percent_mito"]
_ = ax.hexbin(x, y, mincnt=1, xscale = "log")
_ = ax.set_xticks([10, 100, 1000])
_ = ax.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
_ = ax.axvline(500, color="red")
_ = ax.axhline(30, color="red")
_ = plt.xlabel("Number of genes")
_ = plt.ylabel("percent mitochondrial UMIs")
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2-1.png)

Next we log-normalized the data, selected highly variable features and performed PCA

```python
pg.filter_data(all_data)
```

```python
pg.log_norm(all_data)
```

```python
pg.highly_variable_features(all_data, consider_batch = False)
```

```python
pg.pca(all_data)
```

![plot of chunk process_data](figure/unnamed-chunk-2-3.png)

To account for technical variability between samples, we used the Harmony algorithm to align the pricniple component scores.
These adjust PCs were used for downstrseam leiden clustering and UMAP dimensionality reduction.  Once we have the data in UMAP space,
we can look at our QC metrics, along with our basic leiden clustering.

```python
# pg.run_harmony(all_data, n_jobs = 5)
# pg.neighbors(all_data, rep = "pca_harmony")
# pg.diffmap(all_data, rep = "pca_harmony")
# pg.leiden(all_data, rep = "pca_harmony")
# pg.umap(all_data, rep = "pca_harmony")
# pg.write_output("/home/nealpsmith/projects/medoff/data/all_data_harmonized.h5ad")

# Read in the harmonized data
harmonized_data = pg.read_input("/home/nealpsmith/projects/medoff/data/all_data_harmonized.h5ad")
```

```
## 2021-10-01 10:06:57,030 - pegasus - INFO - Time spent on 'read_input' = 5.28s.
```

```python
figure = sc.pl.umap(harmonized_data, color = ["n_genes", "percent_mito", "leiden_labels"],
                    cmap = colormap, return_fig = True, show = False, ncols = 2)
figure.set_size_inches(11, 11)
figure
```

![plot of chunk harmonize](figure/unnamed-chunk-2-5.png)


Now that the data is in UMAP space, we can use cannonical markers to try to define major lineages.

```python
lin_genes = ["EPCAM", "CD8A", "IL7R", "LYZ", "MS4A1", "CPA3", "GNLY"]
figure = sc.pl.umap(harmonized_data, color = lin_genes,
           cmap = colormap, ncols = 3, return_fig = True, show = False,
                    wspace = 0.2, hspace = 0.3)
figure.set_size_inches(10, 7)
```

![plot of chunk marker_genes](figure/unnamed-chunk-2-7.png)
Using all of this info (along with other DEG info), we can assign the clusters to our major lineages.

```python
cell_clust_dict = {
        "1" : "CD8 T cells",
        "2" : "CD4 T cells",
        "3" : "Epithelial",
        "4" : "MPS",
        "5" : "Epithelial",
        "6" : "Epithelial",
        "7" : "Epithelial",
        "8" : "CD8 T and NK cells",
        "9" : "Epithelial",
        "10" : "MPS",
        "11" : "B cells",
        "12" : "CD8 T cells",
        "13" : "MPS",
        "14" : "MPS",
        "15" : "CD8 T cells",
        "16" : "Mast cells",
        "17" : "Epithelial",
        "18" : "Epithelial",
        "19" : "Epithelial",
        "20" : "Epithelial"
    }
harmonized_data.obs["cell_set"] = [cell_clust_dict[clust] for clust in harmonized_data.obs["leiden_labels"]]

figure = sc.pl.umap(harmonized_data, color = "cell_set", return_fig = True, show = False, legend_loc = "on data")
```

```
## ... storing 'cell_set' as categorical
```

```python
figure.set_size_inches(5, 5)
figure
```

![plot of chunk lineage_umap](figure/unnamed-chunk-2-9.png)

Looking at the major lineage markers in a dot plot, we can appreciate how specific they are for our new assignments.

```python
harmonized_data.obs["cell_set"] = pd.Categorical(harmonized_data.obs["cell_set"],
                                                 categories = ["Epithelial", "CD8 T cells", "CD4 T cells",
                                                               "MPS", "B cells", "Mast cells", "CD8 T and NK cells"])

plot = sc.pl.dotplot(harmonized_data, lin_genes, groupby="cell_set",
                     show=False, return_fig=True, title="lineage markers",
                     cmap=colormap, standard_scale = "var")
```

```
## /home/nealpsmith/.conda/envs/sc_analysis/lib/python3.7/site-packages/pandas/core/arrays/categorical.py:2487: FutureWarning: The `inplace` parameter in pandas.Categorical.remove_unused_categories is deprecated and will be removed in a future version.
##   res = method(*args, **kwargs)
```

```python
axes_dict = plot.get_axes()
axes_dict["mainplot_ax"].set_axisbelow(True)
axes_dict["mainplot_ax"].grid()
_ = axes_dict["color_legend_ax"].set_title("scaled expression")
figure = plt.gcf()
figure.set_size_inches(8, 6)
figure
```

![plot of chunk dotplot](figure/unnamed-chunk-2-11.png)

Looking at the embedding density across UMAP space, we can see there are biases in what types of cells were recovered from different sample types.

```python
harmonized_data.obs["pheno_tmpt"] = ["_".join([pheno, tmpt]) for pheno, tmpt in zip(harmonized_data.obs["phenotype"], harmonized_data.obs["sample"])]
sc.tl.embedding_density(harmonized_data, basis = "umap", groupby = "pheno_tmpt", key_added = "pheno_tmpt_dens")
```

```
## ... storing 'pheno_tmpt' as categorical
```

```python
figure = sc.pl.embedding_density(harmonized_data, basis = "umap", key = "pheno_tmpt_dens",
                                 return_fig = True, show = False,
                                 wspace = 0.2, hspace = 0.4, ncols = 3)
figure.set_size_inches(10, 7)
```

![plot of chunk embedding_density](figure/unnamed-chunk-2-13.png)
