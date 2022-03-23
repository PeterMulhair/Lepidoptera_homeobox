import sys
import glob
import json
from ete3 import Tree
from subprocess import call as unix
from joblib import Parallel, delayed

'''
Script to summarize reapeatmasker
output for all TE types for plotting.
 
'''

TE_types = ['DNA transposons', 'Rolling-circles', 'Penelope', 'LTR', 'LINEs', 'SINEs', 'Satellites', 'Simple', 'Low complexity', 'Unclassified']

#Get genome size for species
with open('sp_genome_size.json') as f:
    sp_genome = json.load(f)

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

#Get list of all species names with repeat output
sp_list = []
for sp_TE in glob.glob('~/TE_annotation/repeatMasker/*'):
    if 'tsv' in sp_TE:
        continue
    elif 'py' in sp_TE:
        continue
    else:
        sp_name = sp_TE.split('/')[-1]
        if sp_name == 'PieBrab':
            sp_name = 'PieBras'
    sp_list.append(sp_name)

#Trim species tree to include species with repeat output        
tree = Tree('lepi_busco_tree.nwk', format = 1)

tree.prune(sp_list)
tree.set_outgroup('TinTrin')
tree.write(format=1, outfile="trimmed_sp_tree.nwk")

#Get reverse order of phylogeny for figure
order_sp = []
with open('sp_order.tsv') as f:
    for line in f:
        sp = line.strip()
        order_sp.append(sp)
order_sp = reversed(order_sp)        
        

#Create attributes tsv file with info on TE content
sp_list = []
with open('lep_TE_summ.tsv', 'w') as outF:
    outF.write('Species\tGenome_size\tShx_count\t')
    #for TE_type in sorted(TE_types):
    for TE_type in TE_types:
        outF.write(TE_type + '\t')
    outF.write('\n')

    for sp_name in order_sp:
        if sp_name == 'PieBras':
            sp_name = 'PieBrab'

        outF.write(sp_name + '\t')

        for genome, size in sp_genome.items():
            if sp_name in genome:
                genome_size = size
        outF.write(str(genome_size) + '\t')

        Shx_count = sp_Shx_count[sp_name]
        outF.write(str(Shx_count) + '\t')
        
        for TE_out in glob.glob('../../repeatMasker/' + sp_name + '/*.tbl'):
            TE_perc = {}
            with open(TE_out) as f:
                for line in f:
                    for TE_type in TE_types:
                        if TE_type in line:
                            line = line.strip()
                            lines = line.split(' ')
                            TE_per = lines[-2]
                            TE_per = TE_per.strip('%')
                            TE_perc[TE_type] = TE_per
            for TE_type in TE_types:
                TE_cont = TE_perc[TE_type]
                outF.write(TE_cont + '\t')
        outF.write('\n')
                            


