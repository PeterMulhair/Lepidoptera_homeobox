# Lepidoptera homeobox evolution
This repository contains the datasets and scripts needed to reproduce the results in X

Each folder contains data and code required to recreate results and figures in each section of the manuscript.

If you use any scripts from this repository please cite:
X

## Instructions

* `00_all_homeobox/` contains files on the homeobox genes present in each of the species for each of the homeobox gene classes. It also contains `.tsv` files and `hbx_count_heatmap.R` required to reproduce Figure 1 in the manuscript.
* `01_Hox_gene_cluster/` contains files on genes present in the Hox cluster for each species as well as `plot_Hox_cluster.R` required to reproduce Supplementary Figure 1. It also contains two subdirectories `TAD_analysis/` and `gene_tree/` which contain all necessary data and code required to reproduce analyses and figures for Figure 3 & 4.
* `02_Shx_duplications/` contains files on LINE density in the Hox cluster and its association with Shx gene duplication. All data and code required to reproduce Figure 5 are present. It also contains the subdirectory `TE_annotation/` which contains all code required to annotate TE content in the genomes, as well as to recreate Supplementary figure 3.
* `03_NK_gene_cluster/` contains files on genes present in the NK cluster for each species as well as `plot_NK_cluster.R` required to reproduce Supplementary Figure 2.

---

Genomes used in this analysis from the [Darwin Tree of Life project](https://www.darwintreeoflife.org/) can be downloaded by using code from [here](https://github.com/PeterMulhair/DToL_insects)

Annotation of homeobox genes from all classes using these genomes can be carried out using the [HbxFinder pipeline](https://github.com/PeterMulhair/HbxFinder)

---

<div align="center">
<p align="center">
<img src="https://github.com/PeterMulhair/Lepidoptera_homeobox/blob/main/01_Hox_gene_cluster/figures/Hox_summary_microp.png" width="700" height="350">
</p>
</div>
