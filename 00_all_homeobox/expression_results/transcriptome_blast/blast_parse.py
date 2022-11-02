import json
import argparse
from Bio import SeqIO
import statistics
from collections import defaultdict
from subprocess import call as unix

parse = argparse.ArgumentParser()

parse.add_argument("-s", "--species",type=str, help="name of species to parse",required=True)

args = parse.parse_args()

#Run first blast search
unix(' tblastn -query homeobox.fasta -db ' + args.species + '_transcriptome -outfmt "6 qseqid sseqid stitle evalue pident bitscore qstart qend qlen sstart send slen" -evalue 1e-5 -max_target_seqs 1 -out ' + args.species + '_transcriptome_homeobox.tsv', shell=True)


#Get information on tpm per transcript     
tpm_ID = {}
with open('../kallisto/' + args.species + '_output/abundance.tsv') as f:
    next(f)
    for line in f:
        lines = line.split('\t')
        ID = lines[0]
        tpm = lines[-1].strip()
        tpm = float(tpm)
        tpm_ID[ID] = tpm

#Store info on transcript ID to sequence     
transcript_seq = {}
with open('../trinity_assemblies/trinity_out_dir_' + args.species + '/Trinity.fasta') as f:
    for record in SeqIO.parse(f, 'fasta'):
        ID = record.id
        seq = str(record.seq)
        transcript_seq[ID] = seq

#Print summary info on mean + median tpm values 
median_tmp = statistics.median(tpm_ID.values())
mean_tpm = statistics.mean(tpm_ID.values())
print(median_tmp, mean_tpm)
    
with open('hbx_naming.json') as f:
    hbx_naming = json.load(f)
#print(hbx_naming)

all_hbx_names = {}
hbx_gene_fam = {}
for k, v in hbx_naming.items():
    for gene, names in v.items():
        all_hbx_names[gene] = names
        hbx_gene_fam[names] = k
#print(all_hbx_names)    

hit_genes = defaultdict(list)
with open(args.species + '_transcriptome_homeobox.tsv') as f:
    for line in f:
        lines = line.split('\t')
        query = lines[0]
        hit = lines[1]
        hit_genes[hit].append(query)
        pident = lines[3]
        #if float(pident) < 50:
            #print(line.strip())


with open(args.species + '_homeobox_transcripts.fasta', 'w') as outF:            
    for k, v in hit_genes.items():
        tpm = tpm_ID[k]
        gene_names_list = []
        for gene in v:
            gene_name = gene.split('|')[1]
            try:
                gene_nom = all_hbx_names[gene_name]
            except:
                gene_nom = gene_name
            gene_names_list.append(gene_nom)
        #print(k, v)
        gene_names_list_combined = '_'.join(set(gene_names_list))
        #print(k, set(gene_names_list), tpm)
        trin_seq = transcript_seq[k]
        outF.write('>' + k + '|' + gene_names_list_combined + '\n' + trin_seq + '\n')
    
    
unix('blastx -query ' + args.species + '_homeobox_transcripts.fasta -db homeobox -outfmt "6 qseqid sseqid stitle evalue pident bitscore qstart qend qlen sstart send slen" -evalue 1e-5 -max_target_seqs 1 -out ' + args.species + '_homeobox_recip.tsv', shell=True)

recip_hit_genes = {}
with open(args.species + '_homeobox_recip.tsv') as f:
    for line in f:
        lines = line.split('\t')
        query = lines[0]
        hit = lines[1]
        tr_query = query.split('|')[0]
        recip_hit_genes[tr_query] = hit

tr_name_tpm = defaultdict(list)
for trID, gene in recip_hit_genes.items():
    tpm = tpm_ID[trID]
    gene_names_list = []
    gene_name = gene.split('|')[1]
    try:
        gene_nom = all_hbx_names[gene_name]
    except:
        gene_nom = gene_name

    #print(trID, gene_nom, tpm)

    tr_name_tpm[gene_nom].append(trID.split('_i')[0] + '|' + str(tpm))

with open('homeobox_expression.tsv', 'a+') as outF:
    for hbx in set(all_hbx_names.values()):
        hbx_fam = hbx_gene_fam[hbx]
        #print(hbx)
        if hbx in tr_name_tpm:
            tr_info = tr_name_tpm[hbx]
            tpm_trID = {}
            for info in tr_info:
                tpmval = info.split('|')[1]
                trID = info.split('|')[0]
                tpm_trID[tpmval] = trID
            max_tpm = max(tpm_trID.keys())
            max_trID = tpm_trID[max_tpm]
            #print(hbx, tr_info.split('|')[0], tr_info.split('|')[1])
            #print(hbx, max_trID, max_tpm)
            if float(max_tpm) >= median_tmp:
                outF.write(args.species + '\t' + hbx_fam + '\t' + hbx + '\t' + max_trID + '\t' + str(max_tpm) + '\tsign' + '\n')
            else:
                outF.write(args.species + '\t' + hbx_fam + '\t' + hbx + '\t' + max_trID + '\t' + str(max_tpm) + '\tnot_sign' + '\n')
        else:
            #print(hbx, 'NA', '0')
            outF.write(args.species + '\t' + hbx_fam + '\t' + hbx + '\tNA\t0\tnot_sign\n')
