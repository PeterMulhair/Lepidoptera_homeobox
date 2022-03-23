#!/usr/bin/env python3
import os
import re
import glob
import argparse
import subprocess
from ete3 import Tree
from Bio import SeqIO
from subprocess import call as unix
from joblib import Parallel, delayed

# Author: Peter Mulhair
# Date: 20/01/2022

'''
Script to run HyPhy RELAX model
to test for relaxation or 
intensification of selection on
any of the Lycaenidea branches
in the three Shx genes.
'''

parse = argparse.ArgumentParser()

parse.add_argument("-p", "--path",type=str, help="path to dir with genomes",required=True)
parse.add_argument("-g", "--gene",type=str, help="homeobox gene to analyse",required=True)

args = parse.parse_args()

pwd = os.getcwd()

#Get GCAID to species name
GCA_2_sp = {}
for genome in glob.glob(args.path + '*fasta'):
    GCA = genome.split('/')[-1].split('.')[0]
    with open(genome) as f:
        first_line = f.readline()
        sp = first_line.split(' ')[1] + '_' + first_line.split(' ')[2]
        GCA_2_sp[GCA] = sp

#Store list of species to be used in the analysis        
sp_list = []
with open('../sp_lycan_test.txt') as f:
    for line in f:
        sp = line.strip()
        sp_list.append(sp)

#Create output dirs        
os.makedirs('results', exist_ok=True)
os.makedirs('results/RELAX', exist_ok=True)
os.makedirs('results/RELAX/' + args.gene)
os.makedirs('results/RELAX/' + args.gene + '/prot_align')
os.makedirs('results/RELAX/' + args.gene + '/codon_align')

unix('cp ../lep_' + args.gene + '_AA.fasta results/RELAX/' + args.gene, shell=True)
unix('cp ../lep_' + args.gene + '_nuc.fasta results/RELAX/' + args.gene, shell=True)

os.chdir('results/RELAX/' + args.gene)
for fasta in glob.glob('*_AA.fasta'):
    #Ensure the nucleotide and amino acid sequences match for each gene
    unix('python3 ../../../../translate_nuc.py -i ' + fasta, shell=True)

#Align protein sequences    
unix('mv *_AA.fasta prot_align', shell=True)
os.chdir('prot_align/')
for fas in glob.glob('*fasta'):
    gene = fas.split('.')[0]
    unix('mafft --quiet ' + fas + ' > ' + gene + '_aln.fas', shell=True)
os.chdir(pwd)

#Use pal2nal to get codon alignments
for fasta in glob.glob('results/RELAX/' + args.gene + '/*framenuc*'):
    gene = fasta.split('/')[-1].split('_framenuc')[0]
    prot_aln = 'results/RELAX/' + args.gene + '/prot_align/' + gene + '_AA_aln.fas'
    unix('pal2nal.pl ' + prot_aln + ' ' + fasta + ' >> results/RELAX/' + args.gene + '/codon_align/' + gene + '_aln.fasta -output fasta', shell=True)

#Prune species tree to retain species present in gene fasta file    
subprocess.run("bash -c 'source activate vespasian && vespasian infer-gene-trees --warnings --progress results/RELAX/" + args.gene + "/codon_align/ /home/zoo/zool2500/tree_of_life/raw/updated_DToL_data/phylo/busco/busco_analysis/outgroup_busco/super_tree/astral_tree.nwk -o results/RELAX/" + args.gene + "/gene-trees/ && source deactivate'", shell=True)

codon_aln_list = []
for codon_aln in glob.glob('results/RELAX/' + args.gene + '/codon_align/*fasta'):
    codon_aln_list.append(codon_aln)

#List of species in Lycaenidae    
lycanidae_list = ['Lycaena_phlaeas', 'Celastrina_argiolus', 'Glaucopsyche_alexis', 'Plebejus_argus', 'Cyaniris_semiargus', 'Aricia_agestis', 'Lysandra_bellargus', 'Lysandra_coridon']    

os.chdir('results/RELAX/' + args.gene + '/gene-trees')
with open('../../../../lycan_sp.txt') as f, open('lycan_gene.txt', 'w') as outF:
    for line in f:
        line = line.strip()
        outF.write(line + '|' + args.gene + '\n')

#Root pruned species tree        
t = Tree('lep_' + args.gene + '_aln.nwk')
anc_gene = 'YpsScab|' + args.gene
t.set_outgroup(anc_gene)
t.write(format=1, outfile='lep_' + args.gene + '_aln_rooted.nwk')

#Label all branches related to Lycaenidae inlcuding tip and internal branches
subprocess.run("bash -c 'source activate hyphy && hyphy ../../../../label-tree.bf --tree lep_" + args.gene + "_aln_rooted.nwk --list lycan_gene.txt --output lycan_label.nwk && source deactivate'", shell=True)
os.chdir(pwd)

#Run HyPhy RELAX on coding alignments
def run_hyphy(codon_aln):
    subprocess.run("bash -c 'source activate hyphy && hyphy relax --alignment " + codon_aln + " --tree results/RELAX/" + args.gene + "/gene-trees/lycan_label.nwk --reference-group Foreground && source deactivate'", shell=True)

    
Parallel(n_jobs=1)(delayed(run_hyphy)(aln) for aln in codon_aln_list)
