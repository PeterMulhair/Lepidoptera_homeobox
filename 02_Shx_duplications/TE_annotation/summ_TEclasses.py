#!/usr/bin/env python3
from subprocess import call as unix
import glob
from pathlib import Path
import os

# Author: Peter Mulhair
# Date: 23/08/2021

'''
Short script to parse RepeatMasker
output to pull out TEs from four
major TE classes present on the 
chromosome containing the Hox gene
cluster.
'''

#Major TE classes to parse
TE_classes = ['SINE', 'LINE', 'LTR', 'DNA']

#Store chromosome containing Hox cluster for each species
sp_hoxChr = {}
with open('HOX_cluster_only_filtered.txt') as f:
    lines = f.read()
    sp_info = lines.split('>')
    for info in sp_info[1:]:
        speciesID = info.split('\n')[0]
        spID_short = speciesID.split('_')[2][2:-1]
        if spID_short == 'phHyp':
            spID_short = 'AphHyp'
            
        gene_info = info.split('\n')[1:-1]
        for genes in gene_info:
            if 'Dfd' in genes:
                hbx_chr = genes.split('\t')[0]
                sp_hoxChr[spID_short] = hbx_chr

#Create file of species names and Hox chromosome                
with open('sp_HoxChr.tsv', 'w') as outF:
    for species in glob.glob('*'):
        if '.tsv' in species:
            continue
        elif 'py' in species:
            continue
        else:
            hoxChr = sp_hoxChr[species]
            outF.write(species + '\t' + hoxChr + '\n')

#Parse each species output in repeatmasker output directory            
for sp, chrm in sp_hoxChr.items():
    file_check = Path(sp +'/' + sp + '_HoxChr_DNA.out')
    if file_check.is_file():
        continue
    else:
        try:
            os.chdir(sp)
            print(sp)
            for GCA in glob.glob('GCA*out'):
                GCA_file = GCA

            for TE in TE_classes:
                unix('head -n 2 ' + GCA_file + ' >> ' + sp + '_HoxChr_' + TE + '.out', shell=True)
                unix('grep "' + chrm + '" ' + GCA_file + ' |grep "' + TE + '" >> ' + sp + '_HoxChr_' + TE + '.out', shell=True)

            os.chdir('../')
        except:
            continue
