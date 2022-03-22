library(tidyverse)
library(ggtree)
library(patchwork)
library(phytools)

no_legend <- function() theme(legend.position="none")

tre <- read.newick('lepi_busco_tree_rooted.nwk')
to_drop <- c("LimLuna", "LimMarm", "HypKaha")
p1 <- drop.tip(tre, to_drop)


p <- ggtree(p1, branch.length="none") + geom_tiplab()
p<- ggtree::rotate(p, 124)
p<- ggtree::rotate(p, 142)
p<- ggtree::rotate(p, 143)
p<- ggtree::rotate(p, 178)



#Lepidoptera hbx counts
sp_hbx_count <- read.delim("sp_hbx_count_lepi.tsv")
sp_hbx_genes <- read.delim("sp_hbx_genes_lepi.tsv")


gg_hist <- ggplot(sp_hbx_count, aes(x = factor(Species, levels=unique(Species)), Count)) +
  geom_col(fill = "#446455", alpha=0.6) + no_legend() + coord_flip() + theme_classic() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_text(size = 15)) +
  scale_y_continuous(expand = c(0,0,0.1,0)) + scale_x_discrete(limits = rev)

gg_heat <- ggplot(sp_hbx_genes, aes(fill=Count, x = factor(Gene, levels=unique(Gene)), y = factor(Species, levels=unique(Species)))) + geom_tile() + theme_classic() +
scale_fill_gradientn(breaks= c(1, 2, 3, 4, 5), na.value = '#3B9AB2', lim = c(min(sp_hbx_genes$Count), 5), colours = c('#ffffff', '#fee090', '#E1AF00', '#F21A00', '#78B7C5'), guide = "legend") + theme(axis.text.x = element_text(angle = 70, hjust = 1, size = 6)) + 
theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank()) + scale_y_discrete(limits = rev)

a<-p + gg_heat + gg_hist + plot_annotation(tag_levels="A") & theme(plot.tag = element_text(size = 20))

pdf("hbx_gene_counts_lepi.pdf", height = 15, width = 20)
print(a)
dev.off()


