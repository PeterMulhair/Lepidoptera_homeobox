import toytree
import toyplot
import toyplot.svg
import toyplot.pdf
import numpy as np
import pandas as pd
from ete3 import Tree
from collections import defaultdict
from statsmodels.sandbox.stats.multicomp import multipletests


sp_list = []
with open('sp_list.txt') as f:
    for line in f:
        sp = line.strip()
        sp_list.append(sp)

LINE_sign_sp = []
with open('HOX_LINEdensity_stats.tsv') as f:
    LINE_pvals = []
    species_list = []
    for line in f:
        lines = line.split('\t')
        sp = lines[0]
        pvalue = lines[-1].strip()
        #print(sp, pvalue)
        LINE_pvals.append(float(pvalue))
        species_list.append(sp)
        
    p_adjusted = multipletests(LINE_pvals, method='bonferroni')
    p_adjusted_vals = p_adjusted[1]
    
    sp_count = 0
    for species in species_list:
        p_adjusted_val = p_adjusted_vals[sp_count]
        sp_count+=1
        if float(p_adjusted_val) <= 0.05:
            LINE_sign_sp.append(species)
        else:
            continue

        
#Get Shx loci count for species
sp_Shx_count = {}
with open('HOX_cluster_only_filtered.txt') as f:
    lines = f.read()
    sp_info = lines.split('>')
    for info in sp_info[1:]:
        hox_gene_pos = []
        speciesID = info.split('\n')[0]
        spID_short = speciesID.split('_')[2][2:-1]
        if spID_short == 'phHyp':
            spID_short = 'AphHyp'
            
        gene_info = info.split('\n')[2:-2]#Ignore lab and Ro from Hox cluster
        
        #Get number of Shx genes in cluster
        Shx_count = 0
        shx_gene_pos = []
        for genes in gene_info:
            gene = genes.split('\t')[2]
            if 'Shx' in gene:
                Shx_count += 1
            
        sp_Shx_count[spID_short] = Shx_count

tandDup_sp = []
for sp, shx_count in sp_Shx_count.items():
    if int(shx_count) > 6:
        tandDup_sp.append(sp)
        
#Save info on which species have tandem duplications and enrichment of LINE elements
with open('dup_LINE_summary.csv', 'w') as outF:
    for species in sp_list:
        if species in tandDup_sp:
            outF.write('1,')
        else:
            outF.write('0,')
        if species in LINE_sign_sp:
            outF.write('1\n')
        else:
            outF.write('0\n')

            
sptree = Tree('lepi_busco_tree.nwk', format = 1)

sptree.prune(sp_list)
ancestor = sptree.get_common_ancestor("TinTrin","TinSemi")
sptree.set_outgroup(ancestor)
sptree.write(format=1, outfile="rooted_sp_tree.nwk")

#Hox cluster length data
spdata = np.genfromtxt('HOX_cluster_chr_lengths.tsv', delimiter='\t')

#Hox intergenic dists
#distdata = {}
distdata = defaultdict(list)
with open('HOX_cluster_chr_intergene_dist.tsv') as f:
    next(f)
    for line in f:
        lines = line.split('\t')
        sp = lines[0]
        dist = lines[1].strip()
        distdata[sp].append(dist)
        
tree = toytree.tree("rooted_sp_tree.nwk")

canvas = toyplot.Canvas(width=1000, height=1400)

ax0 = canvas.cartesian(bounds=(80, 200, 50, 1300), padding=15, ymin=0, ymax=122)
ax1 = canvas.cartesian(bounds=(225, 325, 50, 1300), padding=15, ymin=0, ymax=122)

tree.draw(axes=ax0, width=400, height=900, use_edge_lengths=False, tip_labels_style={"-toyplot-anchor-shift": "15px"});
ax0.show = False

spdata = spdata[::-1]
ax1.bars(spdata[:, 1], along='y',color="#fec44f");
ax1.show = True
ax1.y.show = False
ax1.x.ticks.show = True
ax1.x.label.text='Hox cluster size'


colourlist = ['#2c7fb8', '#a63603']
sumdata = np.genfromtxt('dup_LINE_summary.csv', delimiter=',')
ncols = 2
xoffset = 60
for col in range(2):
    data = sumdata[:, col]
    ax0.scatterplot(np.repeat(col, tree.ntips) + xoffset, np.arange(tree.ntips), marker='s', size=10, color=colourlist[col], opacity=0.1 + data[::-1] / data.max(), title=data,);
ax0.x.domain.max = 70
                    
toyplot.svg.render(canvas, 'cluster_stats_leps_info.svg')
