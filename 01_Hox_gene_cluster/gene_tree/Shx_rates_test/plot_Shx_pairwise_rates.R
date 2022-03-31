library(ggplot2)
library(patchwork)
library(getopt)

Shx_pairwiseID <- read.delim("~/Dropbox/Postdoc/tree_of_life/results/gene_cluster_analysis/homeobox_genes/hbxfinder_output/updated_genomes/seq_files/rates_test/Shx_pairwiseID_split.tsv", header=FALSE)

Shx<-ggplot(Shx_pairwiseID, aes(x = V3, y = V1)) + geom_boxplot(aes(fill = V2)) + theme_classic() + scale_fill_manual(values = c("#ff7f00", "#ffdbbb")) +
ylab(NULL) + 
xlab("Percent identity") + 
theme(axis.text.y = element_text(size = 20), axis.title.y = element_blank(), axis.title.x = element_text(size = 20), axis.text.x = element_text(size = 15)) + theme(legend.position = "none")

pdf("lycanidae_pairwise_Shx.pdf", height = 5, width = 7)
print(Shx)
dev.off()




ShxA_pairwiseID <- read.delim("~/Dropbox/Postdoc/tree_of_life/results/gene_cluster_analysis/homeobox_genes/hbxfinder_output/updated_genomes/seq_files/rates_test/ShxA_pairwiseID_split.tsv", header=FALSE)
ShxB_pairwiseID <- read.delim("~/Dropbox/Postdoc/tree_of_life/results/gene_cluster_analysis/homeobox_genes/hbxfinder_output/updated_genomes/seq_files/rates_test/ShxB_pairwiseID_split.tsv", header=FALSE)
ShxC_pairwiseID <- read.delim("~/Dropbox/Postdoc/tree_of_life/results/gene_cluster_analysis/homeobox_genes/hbxfinder_output/updated_genomes/seq_files/rates_test/ShxC_pairwiseID_split.tsv", header=FALSE)

ShxA_split<-split(ShxA_pairwiseID, ShxA_pairwiseID$V1)
ShxA_lyc<-ShxA_split$LycSp
ShxA_nonlyc<-ShxA_split$NonLycSp

ShxA_stats <- wilcox.test(ShxA_lyc$V2, ShxA_nonlyc$V2, paired=FALSE) 
ShxA_pvalue<-ShxA_stats$p.value

ShxB_split<-split(ShxB_pairwiseID, ShxB_pairwiseID$V1)
ShxB_lyc<-ShxB_split$LycSp
ShxB_nonlyc<-ShxB_split$NonLycSp

ShxB_stats <- wilcox.test(ShxB_lyc$V2, ShxB_nonlyc$V2, paired=FALSE) 
ShxB_pvalue<-ShxB_stats$p.value

ShxC_split<-split(ShxC_pairwiseID, ShxC_pairwiseID$V1)
ShxC_lyc<-ShxC_split$LycSp
ShxC_nonlyc<-ShxC_split$NonLycSp

ShxC_stats <- wilcox.test(ShxC_lyc$V2, ShxC_nonlyc$V2, paired=FALSE) 
ShxC_pvalue<-ShxC_stats$p.value

ShxA_pvalue
ShxB_pvalue
ShxC_pvalue