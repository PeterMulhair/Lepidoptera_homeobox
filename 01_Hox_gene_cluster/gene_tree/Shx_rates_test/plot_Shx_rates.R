require(ape)
require(pegas)
library(tidyr)
require(adegenet)

#'#ff7f00', '#ffdbbb'	ShxA
#'#e31a1c', '#e2c3c5'	ShxB
#'#33a02c', '#acc1aa'	ShxC

obj <- fasta2DNAbin("lep_ShxC_aln.fa", chunk=10)


outgroup <- obj[1, ]
n <- nrow(obj)
cc <- combn(2:n, 2)
FUN <- function(x)
    rr.test(obj[x[1], ], obj[x[2], ], outgroup)$Pval
OUT <- apply(cc, 2, FUN)
### two ways to arrange the output:
RES <- matrix(NA, n - 1, n - 1)
RES[row(RES) > col(RES)] <- OUT
RES <- t(RES)
RES[row(RES) > col(RES)] <- OUT
RES <- t(RES)
dimnames(RES) <- list(2:n, 2:n)
RES <- as.dist(RES)

seq_dist<-as.data.frame(as.matrix(RES))



lycan_vals<-seq_dist[c('40', '27', '46', '49', '28', '12', '41', '34')]
lycan_vals<-na.omit(lycan_vals)
nonlycan_vals<-seq_dist[c('2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '29', '30', '31', '32', '33', '35', '36', '37', '38', '39', '42', '43', '44', '45', '47', '48', '50', '51', '52', '53', '54', '55', '56', '57', '58', '59', '60', '61', '62', '63', '64', '65', '66', '67', '68', '69', '70')]
nonlycan_vals<-na.omit(nonlycan_vals)

lycan_data<-lycan_vals %>% gather(species, pvalue)
nonlycan_data<-nonlycan_vals %>% gather(species, pvalue)


all_lep <-rbind(lycan_data, nonlycan_data)
all_lep <- bind_rows(lycan_data, nonlycan_data, .id = 'set')
library(dplyr)

lyc<-ggplot(lycan_data[which(lycan_data$pvalue>0),], aes(x = pvalue, y = species)) + geom_boxplot(aes(fill = '#de2d26')) + theme_classic()
non_lyc<-ggplot(nonlycan_data[which(nonlycan_data$pvalue>0),], aes(x = pvalue, y = species)) + geom_boxplot(aes(fill = '#3182bd')) + theme_classic()



sp_groups <- c("Lycaenidae", "Non-Lycaenidae")

all_lep_plotC<-ggplot(all_lep[which(all_lep$pvalue>0),], aes(fill=set, x = pvalue, y = factor(set, levels=unique(set)))) + 
geom_boxplot() + 
theme_classic() + 
scale_fill_manual(values = c('#33a02c', '#a0bf9d')) + 
ylab("Group") + 
theme(axis.text.y = element_text(size = 20), axis.title.y = element_blank(), axis.title.x = element_text(size = 20), axis.text.x = element_text(size = 15)) +
scale_y_discrete(labels= sp_groups) + theme(legend.position = "none") + expand_limits(x = 0)



pdf("lycanidae_rates_ShxC.pdf", height = 3, width = 15)
print(all_lep_plotC)
dev.off()

