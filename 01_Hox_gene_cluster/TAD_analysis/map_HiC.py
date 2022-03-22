#!/usr/bin/env python3
import os
import glob
import argparse
from Bio import SeqIO
from subprocess import call as unix
from joblib import Parallel, delayed

parse = argparse.ArgumentParser()

parse.add_argument("--path",type=str, help="path to genomes in fasta format",required=True)

args = parse.parse_args()


##Map Hi-C reads to genome 
os.makedirs('../data/results/mapped_files', exist_ok=True)

def map_reads(sp):
    read1 = '../data/raw/' + sp + '_1.fastq.gz'
    read2 = '../data/raw/' + sp + '_2.fastq.gz'

    sp_genome = glob.glob(args.path + '*' + sp + '*')
    sp_genome = sp_genome[0]
    #Make bwa index of genome
    unix('bwa index ' + sp_genome + ' -p ../data/results/' + sp + '_index', shell=True)
    
    #Map Hi-C reads to genome
    unix('bwa mem -t 20 -A1 -B4 -E50 -L0 ../data/results/' + sp + '_index ' + read1 + ' | samtools view -Shb - > ../data/results/mapped_files/' + sp + '_1.bam', shell=True)
    unix('bwa mem -t 20 -A1 -B4 -E50 -L0 ../data/results/' + sp + '_index ' + read2 + ' | samtools view -Shb - > ../data/results/mapped_files/' + sp + '_2.bam', shell=True)

readIDs_sp = {}
sp_list = []

#Path to raw HiC reads in fastq format
for fastq in glob.glob('../data/raw/*fastq.gz'):
    sp = fastq.split('/')[-1].split('_')[0]
    sp_list.append(sp)
    
sp_list = set(sp_list)        

Parallel(n_jobs=17)(delayed(map_reads)(k) for k in sp_list)


##Build HiC matrix and annotate TADs 
os.makedirs('../data/results/TADboundaries', exist_ok=True)
def hic_map(sp):
    unix('hicBuildMatrix --samFiles ../data/results/mapped_files/' + sp + '_1.bam ../data/results/mapped_files/' + sp + '_2.bam --binSize 5000 --inputBufferSize 100000 --threads 10 -o ../data/results/TAD_boundaries/' + sp + '_hic_matrix_5kb.h5 --QCfolder ../data/results/TAD_boundaries/' + sp + '_hic_matrix_5kb_hicQC', shell=True)
      
   
    unix('hicFindTADs -m ../data/results/TAD_boundaries/' + sp + '_hic_matrix_5kb.h5 --correctForMultipleTesting fdr --outPrefix ../data/results/TAD_boundaries/' + sp + '_HiCmatrix_minBouDis20Kb_step2000_thres0.001_delta0.01_resolution.fdr --thresholdComparisons 0.001 --delta 0.01 --minBoundaryDistance 20000 --numberOfProcessors 10', shell=True)
    

Parallel(n_jobs=27)(delayed(hic_map)(k) for k in sp_list)


def plot_hic(sp):
    with open('HOX_cluster_only_filtered_lepi.txt') as f, open('../data/results/TAD_boundaries/' + sp + '_HoxGenes.bed', 'w') as outF:
        lines = f.read()
        sps = lines.split('>')
        for species in sps[1:]:
            sp_line = species.split('\n')
            spID = sp_line[0]
            if sp in spID:
                for gene in sp_line[1:-1]:
                    info = gene.split('\t')
                    chrm = info[0]
                    geneID = info[2]
                    loc = info[1]
                    start = loc.split(', ')[0].strip('(')
                    end = loc.split(', ')[1].strip(')')
                    outF.write(chrm + '\t' + start + '\t' + end + '\t' + geneID + '\n')
                    if geneID == 'Scr':
                        hox_chr = info[0]
    with open('../data/results/TAD_boundaries/hic_track_' + sp + '.ini', 'w') as outF1:
        outF1.write('[x-axis]\nfontsize=10\n\n[hic]\nfile = ' + sp + '_hic_matrix_5kb.h5\ntitle = Threshold 0.05\ncolormap = Spectral_r\ndepth = 18000000\nmin_value = 1\nmax_value = 80\ntransform = log1p\nfile_type = hic_matrix\nshow_masked_bins = false\n\n[tads]\nfile = ' + sp + '_HiCmatrix_minBouDis20Kb_step2000_thres0.001_delta0.01_resolution.fdr_domains.bed\nfile_type = domains\nborder_color = black\ncolor = none\noverlay_previous = share-y\n\n[spacer]\nheight = 0.1\n\n[genes]\nfile = ' + sp + '_HoxGenes.bed\ntitle = Hox genes\ncolor = black\nheight = 18\nlabels = true\nfile_type = bed\n')

    sp_genome = glob.glob(args.path + '*' + sp + '*')
    sp_genome = sp_genome[0]
    with open(sp_genome) as f1:
        for record in SeqIO.parse(f1, 'fasta'):
            ID = record.id
            if ID == hox_chr:
                seq = str(record.seq)
                len_seq = len(seq)
    unix("hicPlotTADs --tracks ../data/results/TAD_boundaries/hic_track_" + sp + ".ini --region " + hox_chr + ":1-" + str(len_seq) + " -t '" + sp + " Hox region TADs' -o ../data/results/TAD_boundaries/" + sp + "_TADs_HoxChrm.pdf", shell=True)

    
Parallel(n_jobs=27)(delayed(plot_hic)(k) for k in sp_list)

print(sp_list)
