Figure 7 : CellPhoneDB
================

``` r
library(reticulate)
use_python("/home/nealpsmith/.conda/envs/sc_analysis/bin/python")
```

In figure 7, we wanted to look at potential cell:cell interactions. We first did this using `CellPhoneDB`

# CellPhoneDB

When using CellPhonedb, we isolated the AA and AC cells at allergen seperately. Then we ran cellphonedb as follows :

<pre>
<p>
#!/bin/bash
cellphonedb method statistical_analysis meta_data.csv count_mtx.h5ad --output-path aa_ag_high_iter --result-precision 20 --iterations 10000 --threads 16
cellphonedb method statistical_analysis meta_data.csv count_mtx.h5ad --output-path ana_ag_high_iter --result-precision 20 --iterations 10000 --threads 16
</p>
</pre>
With the outputs, we next wanted to aggregate the significant interactions and compare what we see with AA and AC. One drawback of cellphonedb is that it does not have great type-I error control. Therefore, we set additional criteria when determining what interactions to consider significant. The criteria were as follows :

-imperical p-value &lt;0.001 -receptor in &gt; 10% of receptor cells -receptor in &gt; 20 cells total -ligand in &gt; 10% of ligand cells -ligand in &gt; 20 total cells

``` r
library(tidyverse)
library(glue)
library(magrittr)

pvals_aa <- read.csv("/home/nealpsmith/projects/medoff/data/cellphonedb/aa/aa_ag/pvalues.txt", sep = "\t")
pvals_ana <- read.csv("/home/nealpsmith/projects/medoff/data/cellphonedb/ana/ana_ag/pvalues.txt", sep = "\t")

means_aa <- read.csv("/home/nealpsmith/projects/medoff/data/cellphonedb/aa/aa_ag/means.txt", sep = "\t")
means_ana <- read.csv("/home/nealpsmith/projects/medoff/data/cellphonedb/ana/ana_ag/means.txt", sep = "\t")

# This has the ranks
sig_means_aa <- read.csv("/home/nealpsmith/projects/medoff/data/cellphonedb/aa/aa_ag/significant_means.txt", sep = "\t")
sig_means_ana <- read.csv("/home/nealpsmith/projects/medoff/data/cellphonedb/ana/ana_ag/significant_means.txt", sep = "\t")

# Need the deconvolution file to get all genes in a complex
deconv_aa <- read.csv("/home/nealpsmith/projects/medoff/data/cellphonedb/aa/aa_ag/deconvoluted.txt", sep = "\t")
deconv_ana <- read.csv("/home/nealpsmith/projects/medoff/data/cellphonedb/ana/ana_ag/deconvoluted.txt", sep = "\t")

# Also want to read in the percent of cells/number of cells expressing every gene
perc_cells_aa <- read.csv("/home/nealpsmith/projects/medoff/data/cellphonedb/perc_of_cells_for_each_gene_AA.csv",
                          row.names = 1)
# Need to fix the T cell names (I'm dumb and inconsistent)
colnames(perc_cells_aa)[grep("t_cell_", colnames(perc_cells_aa))] <- sapply(colnames(perc_cells_aa)[grep("t_cell_", colnames(perc_cells_aa))], function(x) sub("t_", "t", x))

perc_cells_ana <- read.csv("/home/nealpsmith/projects/medoff/data/cellphonedb/perc_of_cells_for_each_gene_ANA.csv",
                          row.names = 1)
colnames(perc_cells_ana)[grep("t_cell_", colnames(perc_cells_ana))] <- sapply(colnames(perc_cells_ana)[grep("t_cell_", colnames(perc_cells_ana))], function(x) sub("t_", "t", x))


n_cells_aa <- read.csv("/home/nealpsmith/projects/medoff/data/cellphonedb/n_cells_for_each_gene_AA.csv", row.names = 1)
colnames(n_cells_aa)[grep("t_cell_", colnames(n_cells_aa))] <- sapply(colnames(n_cells_aa)[grep("t_cell_", colnames(n_cells_aa))], function(x) sub("t_", "t", x))

n_cells_ana <- read.csv("/home/nealpsmith/projects/medoff/data/cellphonedb/n_cells_for_each_gene_ANA.csv", row.names = 1)
colnames(n_cells_ana)[grep("t_cell_", colnames(n_cells_ana))] <- sapply(colnames(n_cells_ana)[grep("t_cell_", colnames(n_cells_ana))], function(x) sub("t_", "t", x))

# Get the columns of interactions
interact_cols <- colnames(pvals_aa)[grepl("\\.", colnames(pvals_aa))]

# Now get all of the unique cellular subsets
cluster_names <- lapply(interact_cols, function(x) strsplit(x, "\\.")[[1]]) %>%
  unlist(use.names = FALSE) %>%
  unique()


perc_cutoff <- 10
ncell_cutoff <- 20

if (! file.exists("/home/nealpsmith/projects/medoff/data/cellphonedb/all_sig_interactions_by_group_strict_cutoff.csv")){
  all_interactions_df <- data.frame("interacting_pair" = character(),
                                  "cluster_pair" = character(),
                                  "aa_perc_a" = numeric(),
                                  "aa_perc_b" = numeric(),
                                  "aa_ncells_a" = numeric(),
                                  "aa_ncells_b" = numeric(),
                                  "aa_rank" = numeric(),
                                  "aa_mean" = numeric(),
                                  "ana_perc_a" = numeric(),
                                  "ana_perc_b" = numeric(),
                                  "ana_ncells_a" = numeric(),
                                  "ana_ncells_b" = numeric(),
                                  "ana_rank" = numeric(),
                                  "ana_mean" = numeric(),
                                  "group" = character())

    for (pair in interact_cols){
    aa_interact <- pvals_aa %>%
      dplyr::select(interacting_pair, pair)
    aa_interact <- aa_interact[aa_interact[[pair]] == 0,]["interacting_pair"]
    if (nrow(aa_interact) > 0){
       aa_interact$cluster_pair <- pair

    # Lets add the other info
    aa_interact$aa_perc_a <- sapply(aa_interact$interacting_pair, function(p){
      exp_cells <- strsplit(pair, ".", fixed = TRUE)[[1]][1]
      gene <- strsplit(p, "_(?!.*_)", perl=TRUE)[[1]][1] # Get the first gene/complex
      if (grepl(" ", gene)){
        # Get the complex members
        genes <- deconv_aa %>%
          dplyr::filter(complex_name == gene) %>%
          .$gene_name
        perc <- max(perc_cells_aa[genes, exp_cells])
        return(perc)
      }  else if (grepl(":", gene)){
        genes <- deconv_aa %>%
          dplyr::filter(complex_name == gene) %>%
          .$gene_name
        perc <- max(perc_cells_aa[genes, exp_cells])
        return(perc)
      } else {
        perc <- perc_cells_aa[gene, exp_cells]
        return(perc)
      }

    })

    aa_interact$aa_perc_b <- sapply(aa_interact$interacting_pair, function(p){
        exp_cells <- strsplit(pair, ".", fixed = TRUE)[[1]][2]
        gene <- strsplit(p, "_(?!.*_)", perl=TRUE)[[1]][2] # Get the second gene/complex
        if (grepl(" ", gene)){
          # Get the complex members
          genes <- deconv_aa %>%
            dplyr::filter(complex_name == gene) %>%
            .$gene_name %>%
            unique()
          perc <- max(perc_cells_aa[genes, exp_cells])
          return(perc)
        }  else if (grepl(":", gene)){
          genes <- deconv_aa %>%
            dplyr::filter(complex_name == gene) %>%
            .$gene_name
          perc <- max(perc_cells_aa[genes, exp_cells])
          return(perc)
        } else {
          perc <- perc_cells_aa[gene, exp_cells]
          return(perc)
        }

      })

     aa_interact$aa_ncells_a <- sapply(aa_interact$interacting_pair, function(p){
      exp_cells <- strsplit(pair, ".", fixed = TRUE)[[1]][1]
      gene <- strsplit(p, "_(?!.*_)", perl=TRUE)[[1]][1] # Get the first gene/complex
      if (grepl(" ", gene)){
        # Get the complex members
        genes <- deconv_aa %>%
          dplyr::filter(complex_name == gene) %>%
          .$gene_name
        perc <- max(n_cells_aa[genes, exp_cells])
        return(perc)
      }  else if (grepl(":", gene)){
        genes <- deconv_aa %>%
          dplyr::filter(complex_name == gene) %>%
          .$gene_name
        perc <- max(n_cells_aa[genes, exp_cells])
        return(perc)
      } else {
        perc <- n_cells_aa[gene, exp_cells]
        return(perc)
      }

    })

    aa_interact$aa_ncells_b <- sapply(aa_interact$interacting_pair, function(p){
      exp_cells <- strsplit(pair, ".", fixed = TRUE)[[1]][2]
      gene <- strsplit(p, "_(?!.*_)", perl=TRUE)[[1]][2] # Get the second gene/complex
      if (grepl(" ", gene)){
        # Get the complex members
        genes <- deconv_aa %>%
          dplyr::filter(complex_name == gene) %>%
          .$gene_name %>%
          unique()
        perc <- max(n_cells_aa[genes, exp_cells])
        return(perc)
      }  else if (grepl(":", gene)){
        genes <- deconv_aa %>%
          dplyr::filter(complex_name == gene) %>%
          .$gene_name
        perc <- max(n_cells_aa[genes, exp_cells])
        return(perc)
      } else {
        perc <- n_cells_aa[gene, exp_cells]
        return(perc)
      }

    })

    # Now want to add the means/rank
    aa_means <- sig_means_aa %>%
      dplyr::select(interacting_pair, rank, pair) %>%
      dplyr::filter(interacting_pair %in% aa_interact$interacting_pair) %>%
     `colnames<-`(c("interacting_pair", "aa_rank", "aa_mean"))
    aa_interact %<>%
      left_join(aa_means, by = "interacting_pair") %>%
      dplyr::filter(aa_perc_a > perc_cutoff, aa_perc_b > perc_cutoff, aa_ncells_a > ncell_cutoff, aa_ncells_b > ncell_cutoff)

    } else {
      aa_interact <- data.frame("interacting_pair" = character(),
                                    "cluster_pair" = character(),
                                    "aa_perc_a" = numeric(),
                                    "aa_perc_b" = numeric(),
                                    "aa_ncells_a" = numeric(),
                                    "aa_ncells_b" = numeric(),
                                    "aa_rank" = numeric(),
                                    "aa_mean" = numeric())
    }

    # Now the ANA info
    ana_interact <- pvals_ana %>%
      dplyr::select(interacting_pair, pair)
    ana_interact <- ana_interact[ana_interact[[pair]] == 0,]["interacting_pair"]

    if(nrow(ana_interact) > 0){
      ana_interact$cluster_pair <- pair

    # Lets add the other info
    ana_interact$ana_perc_a <- sapply(ana_interact$interacting_pair, function(p){
      exp_cells <- strsplit(pair, ".", fixed = TRUE)[[1]][1]
      gene <- strsplit(p, "_(?!.*_)", perl=TRUE)[[1]][1] # Get the first gene/complex
      if (grepl(" ", gene)){
        # Get the complex members
        genes <- deconv_ana %>%
          dplyr::filter(complex_name == gene) %>%
          .$gene_name
        perc <- max(perc_cells_ana[genes, exp_cells])
        return(perc)
      }  else if (grepl(":", gene)){
        genes <- deconv_ana %>%
          dplyr::filter(complex_name == gene) %>%
          .$gene_name
        perc <- max(perc_cells_ana[genes, exp_cells])
        return(perc)
      } else {
        perc <- perc_cells_ana[gene, exp_cells]
        return(perc)
      }

    })

    ana_interact$ana_perc_b <- sapply(ana_interact$interacting_pair, function(p){
        exp_cells <- strsplit(pair, ".", fixed = TRUE)[[1]][2]
        gene <- strsplit(p, "_(?!.*_)", perl=TRUE)[[1]][2] # Get the second gene/complex
        if (grepl(" ", gene)){
          # Get the complex members
          genes <- deconv_ana %>%
            dplyr::filter(complex_name == gene) %>%
            .$gene_name %>%
            unique()
          perc <- max(perc_cells_ana[genes, exp_cells])
          return(perc)
        }  else if (grepl(":", gene)){
          genes <- deconv_ana %>%
            dplyr::filter(complex_name == gene) %>%
            .$gene_name
          perc <- max(perc_cells_ana[genes, exp_cells])
          return(perc)
        } else {
          perc <- perc_cells_ana[gene, exp_cells]
          return(perc)
        }

      })

     ana_interact$ana_ncells_a <- sapply(ana_interact$interacting_pair, function(p){
      exp_cells <- strsplit(pair, ".", fixed = TRUE)[[1]][1]
      gene <- strsplit(p, "_(?!.*_)", perl=TRUE)[[1]][1] # Get the first gene/complex
      if (grepl(" ", gene)){
        # Get the complex members
        genes <- deconv_ana %>%
          dplyr::filter(complex_name == gene) %>%
          .$gene_name
        perc <- max(n_cells_ana[genes, exp_cells])
        return(perc)
      }  else if (grepl(":", gene)){
        genes <- deconv_ana %>%
          dplyr::filter(complex_name == gene) %>%
          .$gene_name
        perc <- max(n_cells_ana[genes, exp_cells])
        return(perc)
      } else {
        perc <- n_cells_ana[gene, exp_cells]
        return(perc)
      }

    })

    ana_interact$ana_ncells_b <- sapply(ana_interact$interacting_pair, function(p){
      exp_cells <- strsplit(pair, ".", fixed = TRUE)[[1]][2]
      gene <- strsplit(p, "_(?!.*_)", perl=TRUE)[[1]][2] # Get the second gene/complex
      if (grepl(" ", gene)){
        # Get the complex members
        genes <- deconv_ana %>%
          dplyr::filter(complex_name == gene) %>%
          .$gene_name %>%
          unique()
        perc <- max(n_cells_ana[genes, exp_cells])
        return(perc)
      }  else if (grepl(":", gene)){
        genes <- deconv_ana %>%
          dplyr::filter(complex_name == gene) %>%
          .$gene_name
        perc <- max(n_cells_ana[genes, exp_cells])
        return(perc)
      } else {
        perc <- n_cells_ana[gene, exp_cells]
        return(perc)
      }

    })

    # Now want to add the means/rank
    ana_means <- sig_means_ana %>%
      dplyr::select(interacting_pair, rank, pair) %>%
      dplyr::filter(interacting_pair %in% ana_interact$interacting_pair) %>%
     `colnames<-`(c("interacting_pair", "ana_rank", "ana_mean"))
    ana_interact %<>%
      left_join(ana_means, by = "interacting_pair") %>%
      dplyr::filter(ana_perc_a > perc_cutoff, ana_perc_b > perc_cutoff, ana_ncells_a > ncell_cutoff, ana_ncells_b > ncell_cutoff)

    } else {
       ana_interact <- data.frame("interacting_pair" = character(),
                                    "cluster_pair" = character(),
                                    "ana_perc_a" = numeric(),
                                    "ana_perc_b" = numeric(),
                                    "ana_ncells_a" = numeric(),
                                    "ana_ncells_b" = numeric(),
                                    "ana_rank" = numeric(),
                                    "ana_mean" = numeric())
    }

    all_interactions <- aa_interact %>%
        dplyr::full_join(ana_interact, by = c("interacting_pair", "cluster_pair")) %>%
        mutate(group = ifelse(!is.na(aa_rank) & !is.na(ana_rank), "both",
                            ifelse(!is.na(aa_rank), "AA", "ANA")))
    # add it to the big df
    all_interactions_df <- rbind(all_interactions_df, all_interactions)

  }

  write.csv(all_interactions_df, "/home/nealpsmith/projects/medoff/data/cellphonedb/all_sig_interactions_by_group_strict_cutoff.csv")

} else {
  all_interactions_df <- read.csv("/home/nealpsmith/projects/medoff/data/cellphonedb/all_sig_interactions_by_group_strict_cutoff.csv", row.names = 1)
}
```

After aggregating the data, we wanted to see which interacting pairs were most specific to each group. Below we show the top 25 pairs for each disease group. We also show this as a circos diagram (figure 7A). This was done with Javascript (D3) and the code for that can be found in the "figure 7a" folder of the repository.

``` r
expand.grid.unique <- function(x, y, include.equals=TRUE)
{
    x <- unique(x)

    y <- unique(y)

    g <- function(i)
    {
        z <- setdiff(y, x[seq_len(i-include.equals)])

        if(length(z)) cbind(x[i], z, deparse.level=0)
    }

    do.call(rbind, lapply(seq_along(x), g))
}

unique_pairs <- expand.grid.unique(cluster_names, cluster_names) %>%
  as.data.frame() %>%
  `colnames<-`(c("clust_1", "clust_2"))

if (! file.exists("/home/nealpsmith/projects/medoff/data/cellphonedb/n_sig_interactions_by_group_clust_pair_one_way.csv")){
  inter_count_df_oneway <- apply(unique_pairs, 1, function(df){
    cl1 = df[["clust_1"]]
    cl2 = df[["clust_2"]]
    pair_data <- inter_count_df %>%
        dplyr::filter(cluster_pair %in% c(glue("{cl1}.{cl2}"), glue("{cl2}.{cl1}")))
      aa_count <- sum(pair_data[pair_data$group == "AA",]$n_interactions)
      ana_count <- sum(pair_data[pair_data$group == "ANA",]$n_interactions)
      both_count <- sum(pair_data[pair_data$group == "both",]$n_interactions)

    df <- data.frame(cluster_pair = rep(glue("{cl1}.{cl2}"), 3),
                 group = c("AA", "ANA", "both"),
                 n_interactions = c(aa_count, ana_count, both_count))
    return(df)
  }) %>%
    do.call(rbind, .)

  inter_count_df_oneway$interaction_type <- sapply(inter_count_df_oneway$cluster_pair, function(x){
    cl1 <- strsplit(x, "\\.")[[1]][1]
    cl2 <- strsplit(x, "\\.")[[1]][2]

    cl1_lin <- strsplit(cl1, "_")[[1]][1]
    cl2_lin <- strsplit(cl2, "_")[[1]][1]

    if (length(unique(c(cl1_lin, cl2_lin))) > 1){
      return("2_lin")
    } else {
      return("1_lin")
    }

  })

  write.csv(inter_count_df_oneway, "/home/nealpsmith/projects/medoff/data/cellphonedb/n_sig_interactions_by_group_clust_pair_one_way.csv")
} else {
  inter_count_df_oneway <- read.csv("/home/nealpsmith/projects/medoff/data/cellphonedb/n_sig_interactions_by_group_clust_pair_one_way.csv",
                                    row.names = 1)
}

# Now get the ones with the largest differential
top_aa_specific_oneway <- inter_count_df_oneway %>%
  dplyr::filter(group != "both", interaction_type == "2_lin") %>%
  dplyr::select(-interaction_type) %>%
  reshape2::dcast(cluster_pair ~ group) %>%
  replace_na(list(AA = 0, ANA = 0)) %>%
  mutate(diff = AA - ANA) %>%
  arrange(by = desc(diff)) %>%
  .[1:25,] %>%
  .$cluster_pair
```

    ## Using n_interactions as value column: use value.var to override.

``` r
# Now get the ones with the largest differential
top_ana_specific_oneway <- inter_count_df_oneway %>%
  dplyr::filter(group != "both", interaction_type == "2_lin") %>%
  dplyr::select(-interaction_type) %>%
  reshape2::dcast(cluster_pair ~ group) %>%
  replace_na(list(AA = 0, ANA = 0)) %>%
  mutate(diff = AA - ANA) %>%
  arrange(by = diff) %>%
  .[1:25,] %>%
  .$cluster_pair
```

    ## Using n_interactions as value column: use value.var to override.

``` r
all_pairs_for_plot <- c(top_aa_specific_oneway, top_ana_specific_oneway)

plot_df = inter_count_df_oneway[inter_count_df_oneway$cluster_pair %in% all_pairs_for_plot,]
plot_df <- plot_df %>%
  dplyr::filter(group != "both")
plot_df %<>%
  dplyr::mutate(value = ifelse(group == "ANA", -n_interactions, n_interactions))

# Need to add annotations
annotations <- c("tcell_7" = "CD8 T (GZMK)",
                 "tcell_3" = "quiesCD8 T ",
                 "tcell_1" = "CD8 T (CLIC3)",
                 "tcell_5" = "Tgd (TRDC)",
                 "tcell_9" = "CD8 T (EGR2)",
                 "tcell_6" = "CD4 Treg (FOXP3)",
                 "tcell_8" = "CD4 Th2 (GATA3)",
                 "tcell_10" = "CD4 ThIFNR (ISG15)",
                 "tcell_2" = "CD4 T (CD40LG)",
                 "tcell_4" = "CD4 Th17 (RORA)",
                 "myeloid_11" = "MC1 (CXCL10)",
                 "myeloid_5" = "MC2 (SPP1)",
                 "myeloid_6" = "MC3 (AREG)",
                 "myeloid_7" = "Mac (FABP4)",
                 "myeloid_13" = "quiesMac",
                 "myeloid_2" = "quiesMC",
                 "myeloid_14" = "Cycling (PCLAF)",
                 "myeloid_1" = "MC4 (CCR2)",
                 "myeloid_8" = "Mac2 (A2M)",
                 "myeloid_4" = "pDC (TCF4)",
                 "myeloid_9" = "migDC (CCR7)",
                 "myeloid_10" = "DC1 (CLEC9A)",
                 "myeloid_3" = "DC2 (CD1C)",
                 "myeloid_12" = "AS DC (AXL)",
                 "epithelial_11" = "Early ciliated",
                 "epithelial_2" = "Ciliated",
                 "epithelial_8" = "Mucous-ciliated",
                 "epithelial_9" = "Hillock",
                 "epithelial_12" = "Deuterosomal",
                 "epithelial_10" = "Cycling basal",
                 "epithelial_5" = "Basal",
                 "epithelial_3" = "Suprabasal",
                 "epithelial_6" = "quiesBasal",
                 "epithelial_13" = "Ionocyte",
                 "epithelial_4" = "Goblet",
                 "epithelial_1" = "quiesGoblet",
                 "epithelial_7" = "Club",
                 "epithelial_14" = "Serous",
                 "NK_cell" = "NK cell",
                 "B_cells" = "B cells",
                 "mast_cells" = "mast cells")

plot_df$annotations <- sapply(as.character(plot_df$cluster_pair), function(x){
  c1 <- annotations[strsplit(x, ".", fixed = TRUE)[[1]][1]]
  c2 <- annotations[strsplit(x, ".", fixed = TRUE)[[1]][2]]
  return(glue("{c1}.{c2}"))
})

annot_order_aa <- plot_df %>%
  dplyr::select(annotations, group, n_interactions) %>%
  group_by(annotations) %>%
  dplyr::mutate(diff = lag(n_interactions) - n_interactions) %>%
  drop_na() %>%
  dplyr::filter(diff > 0) %>%
  dplyr::select(annotations, diff) %>%
  dplyr::arrange(by = diff) %>%
  .$annotations %>%
  as.character()

annot_order_ac <- plot_df %>%
  dplyr::select(annotations, group, n_interactions) %>%
  group_by(annotations) %>%
  dplyr::mutate(diff = lag(n_interactions) - n_interactions) %>%
  drop_na() %>%
  dplyr::filter(diff < 0) %>%
  dplyr::select(annotations, diff) %>%
  dplyr::arrange(by = diff) %>%
  .$annotations %>%
  as.character()

plot_df$annotations <- factor(plot_df$annotations, levels = c(annot_order_ac, annot_order_aa))

# Fix old nomenclature
plot_df$group <- as.character(plot_df$group)
plot_df$group[plot_df$group == "ANA"] <- "AC"

ggplot(plot_df, aes(x = annotations, y = value, group = group, fill = group)) +
  geom_bar(stat = "identity") + coord_flip() +
  theme_classic(base_size = 20) +
  scale_y_continuous(labels = abs) +
  scale_fill_manual(values = c("#FF8000", "#40007F")) +
  ylab("# unique interactions")
```

![](figure_7_cellphonedb_files/figure-markdown_github/side_by_side_bar-1.png)

Next, we wanted to focus in on some interactions of biological interest. We curated a list specific to Th2 interactions. Given Th2 cells are predominantly found in the AA subjects, we displayed the significant ones for the AA only in a dot plot. A large point means it meets our statistical criteria listed at the top. The rank is a value produced by CellPhoneDB and represents how specific an interaction is to given clusters. The more cluster pairs an interaction is found to be significant, the lower the rank. We took the -log10(rank) such that darker points represent interactions that were found in less cluster pairs.

``` r
th2_interact_info <- read.csv("/home/nealpsmith/projects/medoff/data/cellphonedb/th2_interactions_for_dotplot_v3.csv")
t_cell_clusts <- c("tcell_8")
other_clusts <- c("myeloid_1", "myeloid_3", "myeloid_5", "myeloid_6", "myeloid_8", "myeloid_11", "myeloid_12",
                  "epithelial_1", "epithelial_3", "epithelial_4", "epithelial_5", "epithelial_7")


# Lets get the Th2:clust first
th2_first_interactions <- th2_interact_info$interaction[grep("tcell_+", th2_interact_info$pair)]


pval_df <- data.frame()

for (cl1 in t_cell_clusts){
  for(cl2 in other_clusts){
    spec_pvals_aa <- pvals_aa %>%
    dplyr::select(interacting_pair, glue("{cl1}.{cl2}")) %>%
    dplyr::filter(interacting_pair %in% th2_first_interactions) %>%
    `colnames<-`(c("interacting_pair", "p_val")) %>%
    mutate(clust_pair = glue("{cl1}.{cl2}"))

    spec_means_aa <- means_aa %>%
      dplyr::select(interacting_pair, glue("{cl1}.{cl2}")) %>%
      dplyr::filter(interacting_pair %in% th2_first_interactions) %>%
      `colnames<-`(c("interacting_pair", "mean")) %>%
      mutate(clust_pair = glue("{cl1}.{cl2}"))

    # Get the ranks as well
    spec_ranks_aa <- sig_means_aa %>%
      dplyr::select(interacting_pair, rank) %>%
      dplyr::filter(interacting_pair %in% th2_first_interactions) %>%
      `colnames<-`(c("interacting_pair", "rank")) %>%
      mutate(clust_pair = glue("{cl1}.{cl2}"))

    all_info_aa <- list(spec_pvals_aa, spec_means_aa, spec_ranks_aa) %>%
      reduce(left_join, by = c("interacting_pair", "clust_pair")) %>%
      mutate(pheno = "AA")

    spec_pvals_ana <- pvals_ana %>%
    dplyr::select(interacting_pair, glue("{cl1}.{cl2}")) %>%
    dplyr::filter(interacting_pair %in% th2_first_interactions) %>%
    `colnames<-`(c("interacting_pair", "p_val")) %>%
    mutate(clust_pair = glue("{cl1}.{cl2}"))

    spec_means_ana <- means_ana %>%
      dplyr::select(interacting_pair, glue("{cl1}.{cl2}")) %>%
      dplyr::filter(interacting_pair %in% th2_first_interactions) %>%
      `colnames<-`(c("interacting_pair", "mean")) %>%
      mutate(clust_pair = glue("{cl1}.{cl2}"))

    spec_ranks_ana <- sig_means_ana %>%
      dplyr::select(interacting_pair, rank) %>%
      dplyr::filter(interacting_pair %in% th2_first_interactions) %>%
      `colnames<-`(c("interacting_pair", "rank")) %>%
      mutate(clust_pair = glue("{cl1}.{cl2}"))

    all_info_ana <-list(spec_pvals_ana, spec_means_aa, spec_ranks_ana) %>%
      reduce(left_join, by = c("interacting_pair", "clust_pair")) %>%
      mutate(pheno = "ANA")

    all_info <- rbind(all_info_aa, all_info_ana)

    pval_df <- rbind(pval_df, all_info)

  }
}

### Okay now lets get the myl:th2 interactions
th2_second_interactions <- th2_interact_info$interaction[grep("_tcell", th2_interact_info$pair)]

for (cl1 in t_cell_clusts){
  for(cl2 in other_clusts){
    spec_pvals_aa <- pvals_aa %>%
    dplyr::select(interacting_pair, glue("{cl2}.{cl1}")) %>%
    dplyr::filter(interacting_pair %in% th2_second_interactions) %>%
    `colnames<-`(c("interacting_pair", "p_val")) %>%
    mutate(clust_pair = glue("{cl2}.{cl1}"))

    spec_means_aa <- means_aa %>%
      dplyr::select(interacting_pair, glue("{cl2}.{cl1}")) %>%
      dplyr::filter(interacting_pair %in% th2_second_interactions) %>%
      `colnames<-`(c("interacting_pair", "mean")) %>%
      mutate(clust_pair = glue("{cl2}.{cl1}"))

    spec_ranks_aa <- sig_means_aa %>%
      dplyr::select(interacting_pair, rank) %>%
      dplyr::filter(interacting_pair %in% th2_second_interactions) %>%
      `colnames<-`(c("interacting_pair", "rank")) %>%
      mutate(clust_pair = glue("{cl2}.{cl1}"))

    all_info_aa <- list(spec_pvals_aa, spec_means_aa, spec_ranks_aa) %>%
      reduce(left_join, by = c("interacting_pair", "clust_pair")) %>%
      mutate(pheno = "AA")

    spec_pvals_ana <- pvals_ana %>%
    dplyr::select(interacting_pair, glue("{cl2}.{cl1}")) %>%
    dplyr::filter(interacting_pair %in% th2_second_interactions) %>%
    `colnames<-`(c("interacting_pair", "p_val")) %>%
    mutate(clust_pair = glue("{cl2}.{cl1}"))

    spec_means_ana <- means_ana %>%
      dplyr::select(interacting_pair, glue("{cl2}.{cl1}")) %>%
      dplyr::filter(interacting_pair %in% th2_second_interactions) %>%
      `colnames<-`(c("interacting_pair", "mean")) %>%
      mutate(clust_pair = glue("{cl2}.{cl1}"))

    spec_ranks_ana <- sig_means_ana %>%
      dplyr::select(interacting_pair, rank) %>%
      dplyr::filter(interacting_pair %in% th2_second_interactions) %>%
      `colnames<-`(c("interacting_pair", "rank")) %>%
      mutate(clust_pair = glue("{cl2}.{cl1}"))

    all_info_ana <-list(spec_pvals_ana, spec_means_aa, spec_ranks_ana) %>%
      reduce(left_join, by = c("interacting_pair", "clust_pair")) %>%
      mutate(pheno = "ANA")

    all_info <- rbind(all_info_aa, all_info_ana)

    pval_df <- rbind(pval_df, all_info)

  }
}

# Now need to go through the rows and adjust those that we are not considering significant because of additional filters
# (see lines 53-54 for those cutoffs)
# cl_pair <- "tcell_8.myeloid_1"
# inter_pair <- "ICAM2_aLb2 complex"

pval_df$p_val <- apply(pval_df,1, function(df){
  cl_pair <- df[["clust_pair"]]
  inter_pair <- df[["interacting_pair"]]
  info <- all_interactions_df %>%
  dplyr::filter(interacting_pair == inter_pair, cluster_pair == cl_pair)
  if (nrow(info) > 0) {
    if (info$group == "both"){
      pval <- df[["p_val"]]
    } else if (info$group == "AA" & df[["pheno"]] == "AA"){
      pval <- df[["p_val"]]
    } else if (info$group == "AA" & df[["pheno"]] == "ANA"){
      pval <- 1
    } else if (info$group == "ANA" & df[["pheno"]] == "ANA"){
      pval <- df[["p_val"]]
    } else if (info$group == "ANA" & df[["pheno"]] == "AA"){
      pval <- 1
    }
  } else {
    pval <- 1
  }
  return(as.numeric(pval))
})

pval_df$neglogp <- -log10(pval_df$p_val + 0.0001)
pval_df$neglogrank <- -log10(pval_df$rank)


# Lets add the group info
pval_df %<>%
  left_join(dplyr::select(th2_interact_info, interaction, group), by = c("interacting_pair" = "interaction"))

# Okay now the dumb stuff...fix the pairs so they match
pval_df$interacting_pair[grep("[myeloid|epithelial]_[0-9]+.tcell", pval_df$clust_pair)] <-
  sapply(pval_df$interacting_pair[grep("[myeloid|epithelial]_[0-9]+.tcell", pval_df$clust_pair, perl = TRUE)], function(x){
    genes <- strsplit(x, "_")[[1]]
    new_pair <- paste(genes[2], genes[1], sep = "_")
  })

pval_df$clust_pair[grep("[myeloid|epithelial]_[0-9]+.tcell", pval_df$clust_pair)] <-
  sapply(pval_df$clust_pair[grep("[myeloid|epithelial]_[0-9]+.tcell", pval_df$clust_pair)], function(x){
    clusts <- strsplit(x, ".", fixed = TRUE)[[1]]
    new_pair <- paste(clusts[2], clusts[1], sep = ".")
  })

row_order <- c("TNF_TNFRSF1A", "TNFSF14_LTBR", "LTA_LTBR", "CXCR6_CXCL16", "IL13_IL13 receptor", "IL4_IL13 receptor", "IL4_IL4 receptor",
               "CSF1_CSF1R", "CCR4_CCL17", "TGFBR3_TGFB1", "CTLA4_CD86", "CD28_CD86", "ICAM2_aLb2 complex",
               "ICOS_ICOSLG", "PDCD1_CD274", "CCR2_CCL2", "CCR4_CCL22", "CD226_PVR", "TIGIT_PVR", "aEb7 complex_CDH1",
               "IL33 receptor_IL33", "LIF_LIFR")

pval_df$interacting_pair <- factor(pval_df$interacting_pair, levels = rev(row_order))
# Want to just do AA for now
plot_df_aa <- pval_df[pval_df$pheno == "AA",]

# Add annotations
plot_df_aa$annotations <- sapply(as.character(plot_df_aa$clust_pair), function(x){
  c1 <- annotations[strsplit(x, ".", fixed = TRUE)[[1]][1]]
  c2 <- annotations[strsplit(x, ".", fixed = TRUE)[[1]][2]]
  return(glue("{c1}.{c2}"))
})

# Want to order the columns in a specific way
col_order <- c("CD4 Th2 (GATA3).Suprabasal", "CD4 Th2 (GATA3).Basal", "CD4 Th2 (GATA3).Club",
               "CD4 Th2 (GATA3).quiesGoblet", "CD4 Th2 (GATA3).Goblet", "CD4 Th2 (GATA3).MC1 (CXCL10)",
               "CD4 Th2 (GATA3).MC2 (SPP1)", "CD4 Th2 (GATA3).MC3 (AREG)", "CD4 Th2 (GATA3).MC4 (CCR2)",
               "CD4 Th2 (GATA3).DC2 (CD1C)", "CD4 Th2 (GATA3).Mac2 (A2M)", "CD4 Th2 (GATA3).AS DC (AXL)")
plot_df_aa$annotations <- factor(plot_df_aa$annotations, levels = col_order)

ggplot(plot_df_aa, aes(x = annotations, y = interacting_pair, size = neglogp, color = neglogrank)) +
  geom_point() +
  scale_size_continuous(name = "-log10(p-value)") +
  scale_color_continuous(low = "#d3d3d3", high = "#FF8000", limits = c(0, 2.2), name = "-log10(rank)") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
```

![](figure_7_cellphonedb_files/figure-markdown_github/th2_dotplot-1.png)

We can also visualize it where the color represents the mean of the receptor:ligand pair (calculated by cellphonedb)

``` r
ggplot(plot_df_aa, aes(x = annotations, y = interacting_pair, size = neglogp, color = mean)) +
  geom_point() +
  scale_size_continuous(name = "-log10(p-value)") +
  scale_color_continuous(low = "#d3d3d3", high = "#FF8000", limits = c(0, max(pval_df$mean))) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
```

![](figure_7_cellphonedb_files/figure-markdown_github/th2_means-1.png)

Next, we wanted to look at some interactions between basal and myeloid cells. Our prior analyses point at a strong role for basal cells to be playing in the pathology of allergic asthma. Here, we curated a list of interactions that we think could be biologically important.

``` r
epi_myl_interact_info <- read.csv("/home/nealpsmith/projects/medoff/data/cellphonedb/basal_myl_interactions_for_dotplot_v4.csv")
epithelial_clusts <- c("epithelial_5")
myeloid_clusts <- c("myeloid_1", "myeloid_3", "myeloid_5", "myeloid_6", "myeloid_8", "myeloid_11")


epi_first_interactions <- epi_myl_interact_info$interaction[grep("epithelial_+", epi_myl_interact_info$pair)]


pval_df <- data.frame()

for (cl1 in epithelial_clusts){
  for(cl2 in myeloid_clusts){
    spec_pvals_aa <- pvals_aa %>%
    dplyr::select(interacting_pair, glue("{cl1}.{cl2}")) %>%
    dplyr::filter(interacting_pair %in% epi_first_interactions) %>%
    `colnames<-`(c("interacting_pair", "p_val")) %>%
    mutate(clust_pair = glue("{cl1}.{cl2}"))

    spec_means_aa <- means_aa %>%
      dplyr::select(interacting_pair, glue("{cl1}.{cl2}")) %>%
      dplyr::filter(interacting_pair %in% epi_first_interactions) %>%
      `colnames<-`(c("interacting_pair", "mean")) %>%
      mutate(clust_pair = glue("{cl1}.{cl2}"))

    # Get the ranks as well
    spec_ranks_aa <- sig_means_aa %>%
      dplyr::select(interacting_pair, rank) %>%
      dplyr::filter(interacting_pair %in% epi_first_interactions) %>%
      `colnames<-`(c("interacting_pair", "rank")) %>%
      mutate(clust_pair = glue("{cl1}.{cl2}"))

    all_info_aa <- list(spec_pvals_aa, spec_means_aa, spec_ranks_aa) %>%
      reduce(left_join, by = c("interacting_pair", "clust_pair")) %>%
      mutate(pheno = "AA")

    spec_pvals_ana <- pvals_ana %>%
    dplyr::select(interacting_pair, glue("{cl1}.{cl2}")) %>%
    dplyr::filter(interacting_pair %in% epi_first_interactions) %>%
    `colnames<-`(c("interacting_pair", "p_val")) %>%
    mutate(clust_pair = glue("{cl1}.{cl2}"))

    spec_means_ana <- means_ana %>%
      dplyr::select(interacting_pair, glue("{cl1}.{cl2}")) %>%
      dplyr::filter(interacting_pair %in% epi_first_interactions) %>%
      `colnames<-`(c("interacting_pair", "mean")) %>%
      mutate(clust_pair = glue("{cl1}.{cl2}"))

    spec_ranks_ana <- sig_means_ana %>%
      dplyr::select(interacting_pair, rank) %>%
      dplyr::filter(interacting_pair %in% epi_first_interactions) %>%
      `colnames<-`(c("interacting_pair", "rank")) %>%
      mutate(clust_pair = glue("{cl1}.{cl2}"))

    all_info_ana <-list(spec_pvals_ana, spec_means_aa, spec_ranks_ana) %>%
      reduce(left_join, by = c("interacting_pair", "clust_pair")) %>%
      mutate(pheno = "ANA")

    all_info <- rbind(all_info_aa, all_info_ana)

    pval_df <- rbind(pval_df, all_info)

  }
}

### Okay now lets get the myl:th2 interactions
epi_second_interactions <- epi_myl_interact_info$interaction[grep("_epithelial", epi_myl_interact_info$pair)]

for (cl1 in epithelial_clusts){
  for(cl2 in myeloid_clusts){
    spec_pvals_aa <- pvals_aa %>%
    dplyr::select(interacting_pair, glue("{cl2}.{cl1}")) %>%
    dplyr::filter(interacting_pair %in% epi_second_interactions) %>%
    `colnames<-`(c("interacting_pair", "p_val")) %>%
    mutate(clust_pair = glue("{cl2}.{cl1}"))

    spec_means_aa <- means_aa %>%
      dplyr::select(interacting_pair, glue("{cl2}.{cl1}")) %>%
      dplyr::filter(interacting_pair %in% epi_second_interactions) %>%
      `colnames<-`(c("interacting_pair", "mean")) %>%
      mutate(clust_pair = glue("{cl2}.{cl1}"))

    spec_ranks_aa <- sig_means_aa %>%
      dplyr::select(interacting_pair, rank) %>%
      dplyr::filter(interacting_pair %in% epi_second_interactions) %>%
      `colnames<-`(c("interacting_pair", "rank")) %>%
      mutate(clust_pair = glue("{cl2}.{cl1}"))

    all_info_aa <- list(spec_pvals_aa, spec_means_aa, spec_ranks_aa) %>%
      reduce(left_join, by = c("interacting_pair", "clust_pair")) %>%
      mutate(pheno = "AA")

    spec_pvals_ana <- pvals_ana %>%
    dplyr::select(interacting_pair, glue("{cl2}.{cl1}")) %>%
    dplyr::filter(interacting_pair %in% epi_second_interactions) %>%
    `colnames<-`(c("interacting_pair", "p_val")) %>%
    mutate(clust_pair = glue("{cl2}.{cl1}"))

    spec_means_ana <- means_ana %>%
      dplyr::select(interacting_pair, glue("{cl2}.{cl1}")) %>%
      dplyr::filter(interacting_pair %in% epi_second_interactions) %>%
      `colnames<-`(c("interacting_pair", "mean")) %>%
      mutate(clust_pair = glue("{cl2}.{cl1}"))

    spec_ranks_ana <- sig_means_ana %>%
      dplyr::select(interacting_pair, rank) %>%
      dplyr::filter(interacting_pair %in% epi_second_interactions) %>%
      `colnames<-`(c("interacting_pair", "rank")) %>%
      mutate(clust_pair = glue("{cl2}.{cl1}"))

    all_info_ana <-list(spec_pvals_ana, spec_means_aa, spec_ranks_ana) %>%
      reduce(left_join, by = c("interacting_pair", "clust_pair")) %>%
      mutate(pheno = "ANA")

    all_info <- rbind(all_info_aa, all_info_ana)

    pval_df <- rbind(pval_df, all_info)

  }
}

pval_df$p_val <- apply(pval_df,1, function(df){
  cl_pair <- df[["clust_pair"]]
  inter_pair <- df[["interacting_pair"]]
  info <- all_interactions_df %>%
  dplyr::filter(interacting_pair == inter_pair, cluster_pair == cl_pair)
  if (nrow(info) > 0) {
    if (info$group == "both"){
      pval <- df[["p_val"]]
    } else if (info$group == "AA" & df[["pheno"]] == "AA"){
      pval <- df[["p_val"]]
    } else if (info$group == "AA" & df[["pheno"]] == "ANA"){
      pval <- 1
    } else if (info$group == "ANA" & df[["pheno"]] == "ANA"){
      pval <- df[["p_val"]]
    } else if (info$group == "ANA" & df[["pheno"]] == "AA"){
      pval <- 1
    }
  } else {
    pval <- 1
  }
  return(as.numeric(pval))
})

pval_df$neglogp <- -log10(pval_df$p_val + 0.0001)
pval_df$neglogrank <- -log10(pval_df$rank)

# Lets add the group info
pval_df %<>%
  left_join(dplyr::select(epi_myl_interact_info, interaction, group), by = c("interacting_pair" = "interaction"))

# Okay now the dumb stuff...fix the pairs so they match
pval_df$interacting_pair[grep("[myeloid]_[0-9]+.epithelial", pval_df$clust_pair)] <-
  sapply(pval_df$interacting_pair[grep("[myeloid]_[0-9]+.epithelial", pval_df$clust_pair, perl = TRUE)], function(x){
    genes <- strsplit(x, "_")[[1]]
    new_pair <- paste(genes[2], genes[1], sep = "_")
  })

pval_df$clust_pair[grep("[myeloid]_[0-9]+.epithelial", pval_df$clust_pair)] <-
  sapply(pval_df$clust_pair[grep("[myeloid]_[0-9]+.epithelial", pval_df$clust_pair)], function(x){
    clusts <- strsplit(x, ".", fixed = TRUE)[[1]]
    new_pair <- paste(clusts[2], clusts[1], sep = ".")
  })

row_order <- c("F11R_aLb2 complex", "ANXA1_FPR1", "ANXA1_FPR2", "ANXA1_FPR3", "CD55_ADGRE5",
               "TNFSF10_TNFRSF10A", "TNFSF10_TNFRSF10B", "SAA1_FPR2", "IL1 receptor_IL1B",
               "IL1 receptor inhibitor_IL1B", "OSMR_OSM", "WNT5A_FZD1", "WNT5A_FZD2", "FZD6_WNT5A",
               "FZD5_WNT5A", "aVb6 complex_TGFB1", "TGFbeta receptor1_TGFB1","TGFB2_TGFbeta receptor2",
               "EGFR_TGFA", "EGFR_AREG", "EGFR_HBEGF")


pval_df$interacting_pair <- factor(pval_df$interacting_pair, levels = rev(row_order))


# Lets make the plots facetted, and by major cell type
basal_df <- pval_df[grep("epithelial_[35]", pval_df$clust_pair),]

basal_df$annotations <- sapply(as.character(basal_df$clust_pair), function(x){
  c1 <- annotations[strsplit(x, ".", fixed = TRUE)[[1]][1]]
  c2 <- annotations[strsplit(x, ".", fixed = TRUE)[[1]][2]]
  return(glue("{c1}.{c2}"))
})

col_order <- c("Basal.MC1 (CXCL10)", "Basal.MC2 (SPP1)", "Basal.MC3 (AREG)", "Basal.MC4 (CCR2)",
               "Basal.DC2 (CD1C)", "Basal.Mac2 (A2M)")
basal_df$annotations <- factor(basal_df$annotations, levels = col_order)

basal_df$group <- factor(basal_df$group, levels = c("Wnt signaling", "TGFb family", "Leukocyte adhesion", "cytokine",
                                                    "EGFR signaling", "TNFSF10 signaling"))

plot_list <- lapply(c("AA", "ANA"), function(pheno){
  if (pheno == "AA"){
    high_val <- "#FF8000"
  } else {
    high_val <- "#40007F"
  }
  if (pheno == "ANA") {
    plot_pheno <- "AC"
  } else {
    plot_pheno <- "AA"
  }
# Want to make a faceted one
  facet_plot <- ggplot(basal_df[basal_df$pheno == pheno,], aes(x = annotations, y = interacting_pair, size = neglogp,
                                                               color = neglogrank)) +
    geom_point() +
    scale_size_continuous(name = "-log10(p-value)") +
    scale_color_continuous(low = "#d3d3d3", high = high_val, limits = c(0, 2.5), name = "-log10(rank)") +
          ggtitle(glue("Basal interactions : {plot_pheno}")) +
    # facet_wrap(~group + pheno, scales = "free_y", ncol = 1) +
    theme_classic(base_size = 20) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

  # convert ggplot object to grob object
  gp <- ggplotGrob(facet_plot)

  # get gtable columns corresponding to the facets (5 & 9, in this case)
  facet_rows <- gp$layout$t[grepl("panel", gp$layout$name)]

  # get the number of unique x-axis values per facet (1 & 3, in this case)
  y_var <- sapply(ggplot_build(facet_plot)$layout$panel_scales_y,
                  function(l) length(l$range$range))

  # change the relative widths of the facet columns based on
  # how many unique x-axis values are in each facet
  gp$heights[facet_rows] <- gp$heights[facet_rows] * y_var
  return(gp)
})
combined_figure <- ggpubr::ggarrange(plotlist = plot_list, nrow = 1, widths = c(1.7, 1.7))
combined_figure
```

![](figure_7_cellphonedb_files/figure-markdown_github/cellphonedb_basal-1.png)

We were also interested in some of these interactions with the goblet cells. Here, we show similar plots, between goblet cells and myeloid cells.

``` r
epi_myl_interact_info <- read.csv("/home/nealpsmith/projects/medoff/data/cellphonedb/goblet_myl_interactions_for_dotplot.csv")
epithelial_clusts <- c("epithelial_4")
myeloid_clusts <- c("myeloid_1", "myeloid_3", "myeloid_5", "myeloid_6", "myeloid_8", "myeloid_11")


epi_first_interactions <- epi_myl_interact_info$interaction[grep("epithelial_+", epi_myl_interact_info$pair)]


pval_df <- data.frame()

for (cl1 in epithelial_clusts){
  for(cl2 in myeloid_clusts){
    spec_pvals_aa <- pvals_aa %>%
    dplyr::select(interacting_pair, glue("{cl1}.{cl2}")) %>%
    dplyr::filter(interacting_pair %in% epi_first_interactions) %>%
    `colnames<-`(c("interacting_pair", "p_val")) %>%
    mutate(clust_pair = glue("{cl1}.{cl2}"))

    spec_means_aa <- means_aa %>%
      dplyr::select(interacting_pair, glue("{cl1}.{cl2}")) %>%
      dplyr::filter(interacting_pair %in% epi_first_interactions) %>%
      `colnames<-`(c("interacting_pair", "mean")) %>%
      mutate(clust_pair = glue("{cl1}.{cl2}"))

    # Get the ranks as well
    spec_ranks_aa <- sig_means_aa %>%
      dplyr::select(interacting_pair, rank) %>%
      dplyr::filter(interacting_pair %in% epi_first_interactions) %>%
      `colnames<-`(c("interacting_pair", "rank")) %>%
      mutate(clust_pair = glue("{cl1}.{cl2}"))

    all_info_aa <- list(spec_pvals_aa, spec_means_aa, spec_ranks_aa) %>%
      reduce(left_join, by = c("interacting_pair", "clust_pair")) %>%
      mutate(pheno = "AA")

    spec_pvals_ana <- pvals_ana %>%
    dplyr::select(interacting_pair, glue("{cl1}.{cl2}")) %>%
    dplyr::filter(interacting_pair %in% epi_first_interactions) %>%
    `colnames<-`(c("interacting_pair", "p_val")) %>%
    mutate(clust_pair = glue("{cl1}.{cl2}"))

    spec_means_ana <- means_ana %>%
      dplyr::select(interacting_pair, glue("{cl1}.{cl2}")) %>%
      dplyr::filter(interacting_pair %in% epi_first_interactions) %>%
      `colnames<-`(c("interacting_pair", "mean")) %>%
      mutate(clust_pair = glue("{cl1}.{cl2}"))

    spec_ranks_ana <- sig_means_ana %>%
      dplyr::select(interacting_pair, rank) %>%
      dplyr::filter(interacting_pair %in% epi_first_interactions) %>%
      `colnames<-`(c("interacting_pair", "rank")) %>%
      mutate(clust_pair = glue("{cl1}.{cl2}"))

    all_info_ana <-list(spec_pvals_ana, spec_means_aa, spec_ranks_ana) %>%
      reduce(left_join, by = c("interacting_pair", "clust_pair")) %>%
      mutate(pheno = "ANA")

    all_info <- rbind(all_info_aa, all_info_ana)

    pval_df <- rbind(pval_df, all_info)

  }
}

### Okay now lets get the myl:th2 interactions
epi_second_interactions <- epi_myl_interact_info$interaction[grep("_epithelial", epi_myl_interact_info$pair)]

for (cl1 in epithelial_clusts){
  for(cl2 in myeloid_clusts){
    spec_pvals_aa <- pvals_aa %>%
    dplyr::select(interacting_pair, glue("{cl2}.{cl1}")) %>%
    dplyr::filter(interacting_pair %in% epi_second_interactions) %>%
    `colnames<-`(c("interacting_pair", "p_val")) %>%
    mutate(clust_pair = glue("{cl2}.{cl1}"))

    spec_means_aa <- means_aa %>%
      dplyr::select(interacting_pair, glue("{cl2}.{cl1}")) %>%
      dplyr::filter(interacting_pair %in% epi_second_interactions) %>%
      `colnames<-`(c("interacting_pair", "mean")) %>%
      mutate(clust_pair = glue("{cl2}.{cl1}"))

    spec_ranks_aa <- sig_means_aa %>%
      dplyr::select(interacting_pair, rank) %>%
      dplyr::filter(interacting_pair %in% epi_second_interactions) %>%
      `colnames<-`(c("interacting_pair", "rank")) %>%
      mutate(clust_pair = glue("{cl2}.{cl1}"))

    all_info_aa <- list(spec_pvals_aa, spec_means_aa, spec_ranks_aa) %>%
      reduce(left_join, by = c("interacting_pair", "clust_pair")) %>%
      mutate(pheno = "AA")

    spec_pvals_ana <- pvals_ana %>%
    dplyr::select(interacting_pair, glue("{cl2}.{cl1}")) %>%
    dplyr::filter(interacting_pair %in% epi_second_interactions) %>%
    `colnames<-`(c("interacting_pair", "p_val")) %>%
    mutate(clust_pair = glue("{cl2}.{cl1}"))

    spec_means_ana <- means_ana %>%
      dplyr::select(interacting_pair, glue("{cl2}.{cl1}")) %>%
      dplyr::filter(interacting_pair %in% epi_second_interactions) %>%
      `colnames<-`(c("interacting_pair", "mean")) %>%
      mutate(clust_pair = glue("{cl2}.{cl1}"))

    spec_ranks_ana <- sig_means_ana %>%
      dplyr::select(interacting_pair, rank) %>%
      dplyr::filter(interacting_pair %in% epi_second_interactions) %>%
      `colnames<-`(c("interacting_pair", "rank")) %>%
      mutate(clust_pair = glue("{cl2}.{cl1}"))

    all_info_ana <-list(spec_pvals_ana, spec_means_aa, spec_ranks_ana) %>%
      reduce(left_join, by = c("interacting_pair", "clust_pair")) %>%
      mutate(pheno = "ANA")

    all_info <- rbind(all_info_aa, all_info_ana)

    pval_df <- rbind(pval_df, all_info)

  }
}

pval_df$p_val <- apply(pval_df,1, function(df){
  cl_pair <- df[["clust_pair"]]
  inter_pair <- df[["interacting_pair"]]
  info <- all_interactions_df %>%
  dplyr::filter(interacting_pair == inter_pair, cluster_pair == cl_pair)
  if (nrow(info) > 0) {
    if (info$group == "both"){
      pval <- df[["p_val"]]
    } else if (info$group == "AA" & df[["pheno"]] == "AA"){
      pval <- df[["p_val"]]
    } else if (info$group == "AA" & df[["pheno"]] == "ANA"){
      pval <- 1
    } else if (info$group == "ANA" & df[["pheno"]] == "ANA"){
      pval <- df[["p_val"]]
    } else if (info$group == "ANA" & df[["pheno"]] == "AA"){
      pval <- 1
    }
  } else {
    pval <- 1
  }
  return(as.numeric(pval))
})

pval_df$neglogp <- -log10(pval_df$p_val + 0.0001)
pval_df$neglogrank <- -log10(pval_df$rank)

# Lets add the group info
pval_df %<>%
  left_join(dplyr::select(epi_myl_interact_info, interaction, group), by = c("interacting_pair" = "interaction"))

# Okay now the dumb stuff...fix the pairs so they match
pval_df$interacting_pair[grep("[myeloid]_[0-9]+.epithelial", pval_df$clust_pair)] <-
  sapply(pval_df$interacting_pair[grep("[myeloid]_[0-9]+.epithelial", pval_df$clust_pair, perl = TRUE)], function(x){
    genes <- strsplit(x, "_")[[1]]
    new_pair <- paste(genes[2], genes[1], sep = "_")
  })

pval_df$clust_pair[grep("[myeloid]_[0-9]+.epithelial", pval_df$clust_pair)] <-
  sapply(pval_df$clust_pair[grep("[myeloid]_[0-9]+.epithelial", pval_df$clust_pair)], function(x){
    clusts <- strsplit(x, ".", fixed = TRUE)[[1]]
    new_pair <- paste(clusts[2], clusts[1], sep = ".")
  })


row_order <- c("F11R_aLb2 complex", "ANXA1_FPR1", "ANXA1_FPR2", "ANXA1_FPR3", "CD55_ADGRE5",
               "TNFSF10_TNFRSF10A", "TNFSF10_TNFRSF10B", "SAA1_FPR2", "IL1 receptor_IL1B",
               "IL1 receptor inhibitor_IL1B", "OSMR_OSM", "WNT5A_FZD1", "WNT5A_FZD2", "FZD6_WNT5A",
               "FZD5_WNT5A", "aVb6 complex_TGFB1", "TGFbeta receptor1_TGFB1","TGFB2_TGFbeta receptor2",
               "EGFR_TGFA", "EGFR_AREG", "EGFR_HBEGF", "VEGFA_NRP2")


pval_df$interacting_pair <- factor(pval_df$interacting_pair, levels = rev(row_order))

# Lets make the plots facetted, and by major cell type
goblet_df <- pval_df[grep("epithelial_[14]", pval_df$clust_pair),]

goblet_df$annotations <- sapply(as.character(goblet_df$clust_pair), function(x){
  c1 <- annotations[strsplit(x, ".", fixed = TRUE)[[1]][1]]
  c2 <- annotations[strsplit(x, ".", fixed = TRUE)[[1]][2]]
  return(glue("{c1}.{c2}"))
})

# Want to order the columns in a specific way
col_order <- c("Goblet.MC1 (CXCL10)", "Goblet.MC2 (SPP1)", "Goblet.MC3 (AREG)",
               "Goblet.MC4 (CCR2)", "Goblet.DC2 (CD1C)", "Goblet.Mac2 (A2M)")
goblet_df$annotations <- factor(goblet_df$annotations, levels = col_order)

goblet_df$group <- factor(goblet_df$group, levels = c("Wnt signaling", "TGFb family", "Leukocyte adhesion", "cytokine", "EGFR signaling", "TNFSF10 signaling"))

plot_list <- lapply(c("AA", "ANA"), function(pheno){
  if (pheno == "AA"){
    high_val <- "#FF8000"
  } else {
    high_val <- "#40007F"
  }

  if (pheno == "ANA") {
    plot_pheno <- "AC"
  } else {
    plot_pheno <- "AA"
  }

# Want to make a faceted one
  facet_plot <- ggplot(goblet_df[goblet_df$pheno == pheno,], aes(x = annotations, y = interacting_pair, size = neglogp,
                                                                 color = neglogrank)) +
    geom_point() +
    scale_size_continuous(name = "-log10(p-value)") +
    scale_color_continuous(low = "#d3d3d3", high = high_val, limits = c(0, 2.5), name = "-log10(rank)") +
          ggtitle(glue("Goblet interactions : {plot_pheno}")) +
    theme_classic(base_size = 20) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

  # convert ggplot object to grob object
  gp <- ggplotGrob(facet_plot)

  # get gtable columns corresponding to the facets (5 & 9, in this case)
  facet_rows <- gp$layout$t[grepl("panel", gp$layout$name)]

  # get the number of unique x-axis values per facet (1 & 3, in this case)
  y_var <- sapply(ggplot_build(facet_plot)$layout$panel_scales_y,
                  function(l) length(l$range$range))

  # change the relative widths of the facet columns based on
  # how many unique x-axis values are in each facet
  gp$heights[facet_rows] <- gp$heights[facet_rows] * y_var
  return(gp)
})
combined_figure <- ggpubr::ggarrange(plotlist = plot_list, nrow = 1, widths = c(1.7, 1.7))
combined_figure
```

![](figure_7_cellphonedb_files/figure-markdown_github/cellphonedb_goblet-1.png)

Again, we can show this where the point color is represented by the mean of the pair (calculated by cellphonedb)

``` r
plot_list <- lapply(c("AA", "ANA"), function(pheno){
  if (pheno == "AA"){
    high_val <- "#FF8000"
  } else {
    high_val <- "#40007F"
  }

  if (pheno == "ANA") {
    plot_pheno <- "AC"
  } else {
    plot_pheno <- "AA"
  }
  # Want to make a faceted one
  facet_plot <- ggplot(goblet_df[goblet_df$pheno == pheno,], aes(x = annotations, y = interacting_pair, size = neglogp,
                                                                 color = mean)) +
    geom_point() +
    scale_size_continuous(name = "-log10(p-value)") +
    scale_color_continuous(low = "#d3d3d3", high = high_val, limits = c(0, 4), name = "mean") +
                    ggtitle(glue("Goblet interactions : {plot_pheno}")) +
    theme_classic(base_size = 20) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

  # convert ggplot object to grob object
  gp <- ggplotGrob(facet_plot)

  # get gtable columns corresponding to the facets (5 & 9, in this case)
  facet_rows <- gp$layout$t[grepl("panel", gp$layout$name)]

  # get the number of unique x-axis values per facet (1 & 3, in this case)
  y_var <- sapply(ggplot_build(facet_plot)$layout$panel_scales_y,
                  function(l) length(l$range$range))

  # change the relative widths of the facet columns based on
  # how many unique x-axis values are in each facet
  gp$heights[facet_rows] <- gp$heights[facet_rows] * y_var
  return(gp)
})
combined_figure <- ggpubr::ggarrange(plotlist = plot_list, nrow = 1, widths = c(1.7, 1.7))
combined_figure
```

![](figure_7_cellphonedb_files/figure-markdown_github/goblet_means-1.png)
