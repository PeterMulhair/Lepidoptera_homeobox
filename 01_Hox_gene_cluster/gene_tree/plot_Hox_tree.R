library(plyr)
library(ggtree)
library(getopt)

tree <- read.tree("lep_HOX_AA.aln.contree")

hbx_tree <- groupClade(tree, c(1251, 1210, 1158, 1347, 1014, 1012))

hbx_geneNames <- read.delim("lep_HOX_geneNames.tsv", header=FALSE)


p <- ggtree(hbx_tree, aes(color=group), layout='circular', size = 0.0001) +
    xlim(-2, NA) +
    theme(legend.position="none") + scale_color_manual(values=c("black", "#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#dfc27d", "#dfc27d"))

p<- ggtree::rotate(p, 1151)


p3 <- open_tree(p, 100)
p4<-rotate_tree(p3, 270)


pdf("HOX_tree_lepi.pdf", height = 10, width = 10)
print(p4)
dev.off()

