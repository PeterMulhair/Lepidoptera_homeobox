library(ggtree)
library(phytools)
library(gggenes)
library(ggplot2)


lep_full <- read.newick("lepi_busco_tree_rooted.nwk")
lept<- compute.brlen(lep_full, 1)

lept<-drop.tip(lept, 'HypKaha')

tree1 <- reroot(lept, node = 214)


p <- ggtree(tree1, branch.length="none") + geom_tiplab()

p<- ggtree::rotate(p, 126)
p<- ggtree::rotate(p, 133)
p<- ggtree::rotate(p, 134)
p<- ggtree::rotate(p, 202)
p<- ggtree::rotate(p, 238)
p<- ggtree::rotate(p, 239)


HOX_cluster_chr_only_filtered_short_lepi <- read.delim("HOX_cluster_only_filtered_short_lepi.tsv")

HOX_cluster_chr_only_filtered_short_lepi$Molecule <- factor(HOX_cluster_chr_only_filtered_short_lepi$Molecule, levels = unique(HOX_cluster_chr_only_filtered_short_lepi$Molecule))
HOX_cluster_chr_only_filtered_short_lepi$gene <- factor(HOX_cluster_chr_only_filtered_short_lepi$gene, levels = unique(HOX_cluster_chr_only_filtered_short_lepi$gene))

HOX_short <- ggplot(HOX_cluster_chr_only_filtered_short_lepi, aes(xmin = start, xmax = end, y = Molecule, fill = gene, label = gene, order= gene)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(0, "mm")) +
  geom_gene_label(align = "left") +
  facet_wrap(~ Molecule, scales = "free", ncol = 1) +
  scale_fill_manual(values = c("lab" = "#65919b", "Abd-B" = "#dfc27d", "abd-A" = "#dfc27d", "Ubx" = "#dfc27d", "Antp" = "#dfc27d", "ftz" = "#dfc27d", "Scr" = "#dfc27d", "Dfd" = "#dfc27d", "zen" = "#dfc27d", "Pb" = "#dfc27d", "ShxD" = "#1f78b4", "ShxC" = "#33a02c", "ShxB" = "#e31a1c", "ShxA" = "#ff7f00", "Ro" = "#6a3d9a")) +
   theme_genes() +
   theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_text(size = 10))


HOX_tre <- p + HOX_short

pdf("HOX_cluster_plot.pdf", height = 30, width = 20)
print(HOX_tre)
dev.off()
