library(ggplot2)
library(dplyr)


homeobox_expression <- read.delim("homeobox_expression.tsv", header=FALSE)
homeobox_expression[homeobox_expression== 0] <- NA  

gene_list<-c('Abd-B', 'abd-A', 'Ubx', 'Antp', 'Ftz', 'Scr', 'Dfd', 'zen', 'ShxA', 'ShxB', 'ShxC', 'ShxD', 'Pb', 'lab', 'Ro', 'eve', 'exex', 'unpg', 'ind', 'btn', 'cad', 'Lbx', 'NK3', 'NK4', 'Msx', 'Hmx', 'Tlx', 'NK1', 'Emx', 'Abox', 'Barhl', 'Bari', 'En', 'Nedx', 'Dlx', 'NK6', 'NK7', 'Hhex', 'NK2.1', 'Bsx', 'Dbx', 'Msxlx', 'Hlx', 'Noto', 'NK2.2', 'Uncx', 'Pitx', 'Gsc', 'Phox', 'Pax4/6', 'Pax2/5/8', 'Hbn', 'Rax', 'Otp', 'CG11294', 'Shox', 'Arx', 'Repo', 'Prrx', 'Vsx', 'Prop', 'Drgx', 'Pax3/7', 'Otx', 'Lhx2/9', 'Lhx6/8', 'Lmx', 'Lhx1/5', 'Lhx3/4', 'Isl', 'Pou3', 'Pou4', 'Pou2', 'Pou6', 'Six3/6', 'Six1/2', 'Six4/5', 'Meis', 'Irx', 'Pbx', 'Tgif', 'Mkx', 'Cmp', 'Onecut', 'Cux', 'Prox', 'Zfhx', 'Cers')

p<-ggplot(homeobox_expression, aes(fill=V5, x = factor(V3, level=gene_list), y = V1, colour = V6)) + geom_tile() + theme_classic() +
scale_colour_manual(values = c("#ffffff", "#990000")) +
scale_fill_gradientn(na.value = '#f0f0f0', colours = c('#d0d1e6', '#3690c0', '#016c59', '#014636'), guide = "legend") + theme(axis.text.x = element_text(angle = 70, hjust = 1, size = 6)) + 
theme(axis.text.y = element_text(size = 10), axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), strip.background = element_blank(), strip.text.y = element_blank(), legend.position = "none", panel.spacing = unit(0.005, "lines")) + scale_y_discrete(limits = rev)

pdf("homeobox_expression.pdf", height = 6, width = 18)
print(p)
dev.off()


######################################################
#Add dots for genes expressed over a certain threshold
######################################################

homeobox_expression <- read.delim("homeobox_expression.tsv", header=FALSE)
homeobox_expression[homeobox_expression== 0] <- NA  

gene_list<-c('Abd-B', 'abd-A', 'Ubx', 'Antp', 'Ftz', 'Scr', 'Dfd', 'zen', 'ShxA', 'ShxB', 'ShxC', 'ShxD', 'Pb', 'lab', 'Ro', 'eve', 'exex', 'unpg', 'ind', 'btn', 'cad', 'Lbx', 'NK3', 'NK4', 'Msx', 'Hmx', 'Tlx', 'NK1', 'Emx', 'Abox', 'Barhl', 'Bari', 'En', 'Nedx', 'Dlx', 'NK6', 'NK7', 'Hhex', 'NK2.1', 'Bsx', 'Dbx', 'Msxlx', 'Hlx', 'Noto', 'NK2.2', 'Uncx', 'Pitx', 'Gsc', 'Phox', 'Pax4/6', 'Pax2/5/8', 'Hbn', 'Rax', 'Otp', 'CG11294', 'Shox', 'Arx', 'Repo', 'Prrx', 'Vsx', 'Prop', 'Drgx', 'Pax3/7', 'Otx', 'Lhx2/9', 'Lhx6/8', 'Lmx', 'Lhx1/5', 'Lhx3/4', 'Isl', 'Pou3', 'Pou4', 'Pou2', 'Pou6', 'Six3/6', 'Six1/2', 'Six4/5', 'Meis', 'Irx', 'Pbx', 'Tgif', 'Mkx', 'Cmp', 'Onecut', 'Cux', 'Prox', 'Zfhx', 'Cers')

p<-ggplot(homeobox_expression, aes(fill=V5, x = factor(V3, level=gene_list), y = V1)) + geom_tile(colour = "#ffffff", size = 0.5) + theme_classic() +
geom_point(aes(size=V6)) +
scale_size_manual(values=c(sign=1, not_sign=NA), guide="none") +
scale_fill_gradientn(na.value = '#f0f0f0', colours = c('#d0d1e6', '#3690c0', '#016c59', '#014636'), guide = "legend") + theme(axis.text.x = element_text(angle = 70, hjust = 1, size = 6)) + 
theme(axis.text.y = element_text(size = 10), axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), strip.background = element_blank(), strip.text.y = element_blank(), legend.position = "none", panel.spacing = unit(0.005, "lines")) + scale_y_discrete(limits = rev)

pdf("homeobox_expression1.pdf", height = 6, width = 12)
print(p)
dev.off()