#!/usr/bin/env python3

import glob
import argparse
from Bio import SeqIO
from collections import defaultdict, Counter
from subprocess import call as unix

parse = argparse.ArgumentParser()

parse.add_argument("-i", "--input",type=str, help="name of input amino acid fasta file",required=True)

args = parse.parse_args()


def translate_codon(dna):
    table = {'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',  'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',  'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',  'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',  'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',  'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',  'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',  'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',  'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',  'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',  'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',  'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',  'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',  'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',  'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',  'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}
    protein = []
    end = len(dna) - (len(dna) %3) - 1
    for i in range(0,end,3):
        codon = dna[i:i+3]
        if codon in table:
            aminoacid = table[codon]
            protein.append(aminoacid)
        else:
            protein.append("N")
    return "".join(protein)


def rev_seq(seq):
    trans=[]
    for i in seq:
        if i=='A':
            trans.append('T')
        elif i=='C':
            trans.append('G')
        elif i=='G':
            trans.append('C')
        elif i=='T':
            trans.append('A')
        else:
            trans.append(i)
    trans=''.join(trans)
    seq_rev= trans[::-1]
    return seq_rev

def six_frames(seq):
    frames = {'+1':[],'+2':[],'+3':[],'-1':[],'-2':[],'-3':[]}
    seq_rev = rev_seq(seq)
    for j in range(0,3):
        temp = ''.join([seq[j::]])
        temp_rev = ''.join([seq_rev[j::]])
        seq_trans = translate_codon(temp)
        seq_rev_trans = translate_codon(temp_rev)
        if j==0:
            frames['+1']=seq_trans
            frames['-1']=seq_rev_trans
        if j==1:
            frames['+2']=seq_trans
            frames['-2']=seq_rev_trans
        if j==2:
            frames['+3']=seq_trans
            frames['-3']=seq_rev_trans
            
    return frames


nuc_file = args.input.split('_AA')[0] + '_nuc.fasta'

sp_AA_seq = {}
with open(args.input) as f:
    for record in SeqIO.parse(f, 'fasta'):
        header = record.description
        seq = str(record.seq)
        sp_AA_seq[header] = seq

issue_sp = []
with open(nuc_file) as f, open(nuc_file.split('_nuc')[0] + '_framenuc.fasta' ,'w') as outF:
    for record in SeqIO.parse(f, 'fasta'):
        header = record.description
        dnaseq = str(record.seq)
        dnaseq = dnaseq.upper()
        AAseq = sp_AA_seq[header]
        
        translate_seq = six_frames(dnaseq)
        for orient, seq in translate_seq.items():
            if (seq in AAseq) or (AAseq in seq):
                frame = orient[1]
                frame = int(frame) -1
                direc = orient[0]
                if direc == '+':
                    temp_seq = ''.join(dnaseq[frame::])
                else:
                    seq_rev = rev_seq(dnaseq)
                    temp_seq = ''.join([seq_rev[frame::]])

                outF.write('>' + header + '\n' + temp_seq + '\n')
                
            else:
                issue_sp.append(header)

issue_sp_count = Counter(issue_sp)
#print(issue_sp_count)
for k, v in issue_sp_count.items():
    if v > 5:
        print(k)
