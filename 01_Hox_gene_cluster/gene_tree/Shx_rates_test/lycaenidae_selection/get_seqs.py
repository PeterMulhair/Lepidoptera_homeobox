import glob
import argparse
from Bio import SeqIO

parse = argparse.ArgumentParser()

parse.add_argument("-p", "--path",type=str, help="path to dir with genomes",required=True)

args = parse.parse_args()

GCA_2_sp = {}
for genome in glob.glob(args.path + '*fasta'):
    GCA = genome.split('/')[-1].split('.')[0]
    with open(genome) as f:
        first_line = f.readline()
        sp = first_line.split(' ')[1] + '_' + first_line.split(' ')[2]
        GCA_2_sp[GCA] = sp

#Open list of species with normal Shx gene count & lycanidae species (missing ShxD)        
sp_list = []
with open('sp_lycan_test.txt') as f:
    for line in f:
        sp = line.strip()
        sp_list.append(sp)


with open('../../HOX_HD_AA.fasta') as f, open('lep_ShxA_AA.fasta', 'w') as outF, open('lep_ShxB_AA.fasta', 'w') as outF1, open('lep_ShxC_AA.fasta', 'w') as outF2:
    for record in SeqIO.parse(f, 'fasta'):
        ID = record.id
        seq = str(record.seq)

        GCAID = ID.split('_')[0] + '_' + ID.split('_')[1]
        spID = GCA_2_sp[GCAID]
        gene = ID.split('_')[2]
        sp_short = spID.split('_')[0][:3] + spID.split('_')[1][:4].title()
        if spID in sp_list:
            if 'ShxA' in gene:
                outF.write('>' + sp_short + '|' + gene + '\n' + seq + '\n')
            elif 'ShxB' in gene:
                outF1.write('>' + sp_short + '|' + gene + '\n' + seq + '\n')
            elif 'ShxC' in gene:
                outF2.write('>' + sp_short + '|' + gene + '\n' + seq + '\n')


with open('../../HOX_HD_nuc.fasta') as f, open('lep_ShxA_nuc.fasta', 'w') as outF, open('lep_ShxB_nuc.fasta', 'w') as outF1, open('lep_ShxC_nuc.fasta', 'w') as outF2:
    for record in SeqIO.parse(f, 'fasta'):
        ID = record.id
        seq = str(record.seq)
        seq = seq.upper()
        
        GCAID = ID.split('_')[0] + '_' + ID.split('_')[1]
        spID = GCA_2_sp[GCAID]
        gene = ID.split('_')[2]
        sp_short = spID.split('_')[0][:3] + spID.split('_')[1][:4].title()
        if spID in sp_list:
            if 'ShxA' in gene:
                outF.write('>' + sp_short + '|' + gene + '\n' + seq + '\n')
            elif 'ShxB' in gene:
                outF1.write('>' + sp_short + '|' + gene + '\n' + seq + '\n')
            elif 'ShxC' in gene:
                outF2.write('>' + sp_short + '|' + gene + '\n' + seq + '\n')
                
