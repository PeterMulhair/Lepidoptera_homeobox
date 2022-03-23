library(phytools)
library(ggtree)
library(gggenes)
library(ggplot2)

lep_full <- read.newick("lepi_busco_tree_rooted.nwk")
lept<- compute.brlen(lep_full, 1)

lept<-drop.tip(lept, 'HypKaha')

tree1 <- reroot(lept, node = 214)


p <- ggtree(tree1, branch.length="none") + geom_tiplab()
#p + geom_nodelab(aes(label = node))

p<- ggtree::rotate(p, 126)
p<- ggtree::rotate(p, 133)
p<- ggtree::rotate(p, 134)
p<- ggtree::rotate(p, 202)
p<- ggtree::rotate(p, 238)
p<- ggtree::rotate(p, 239)



NK_cluster_chr_only_filtered_short_lepi <- read.delim("NK_cluster_only_filtered_short_lepi.tsv")

NK_cluster_chr_only_filtered_short_lepi$Molecule <- factor(NK_cluster_chr_only_filtered_short_lepi$Molecule, levels = unique(NK_cluster_chr_only_filtered_short_lepi$Molecule))
NK_cluster_chr_only_filtered_short_lepi$gene <- factor(NK_cluster_chr_only_filtered_short_lepi$gene, levels = unique(NK_cluster_chr_only_filtered_short_lepi$gene))

NK_short <- ggplot(NK_cluster_chr_only_filtered_short_lepi, aes(xmin = start, xmax = end, y = Molecule, fill = gene, label = gene, order= gene)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(0, "mm")) +
  geom_gene_label(align = "left") +
  facet_wrap(~ Molecule, scales = "free", ncol = 1) +
  scale_fill_manual(values = c("NK6" = "#dfc27d", "NK7" = "#dfc27d", "Tlx" = "#a6cee3", "Msx" = "#1f78b4", "NK4" = "#b2df8a", "NK3" = "#33a02c", "Lbx" = "#fb9a99", "NK1" = "#e31a1c", "Hmx" = "#fdbf6f", "Emx" = "#ff7f00", "Abox" = "#cab2d6", "Bari" = "#6a3d9a", "Prrx" = "#ffff99", "NK2" = "#b15928")) +
   theme_genes() +
   theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_text(size = 10))


NK_tre <- p + NK_short

pdf("NK_cluster_plot.pdf", height = 30, width = 20)
print(NK_tre)
dev.off()
