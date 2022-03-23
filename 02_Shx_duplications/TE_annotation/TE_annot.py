#!/usr/bin/env python3
import os
import glob
from subprocess import call as unix
from joblib import Parallel, delayed

# Author: Peter Mulhair
# Date: 23/08/2021
# Usage python3 TE_annot.py


os.mkdir('CustomLib')
os.mkdir('RepeatModeler')
os.mkdir('protein_coding_filter')
os.mkdir('repeatMasker')

pwd = os.getcwd()

def TE_analysis(genome):
	genome_short = genome.split('.')[1].split('_')[1][2:-1]

	#Run RepeatModeler
	os.chdir('RepeatModeler')
	os.mkdir(genome_short)
	os.chdir(genome_short)

	unix('~/bin/RepeatModeler-2.0.2a/BuildDatabase -name ' + genome_short + ' -engine ncbi ~/TE_annotation/data/genomes/' + genome, shell=True)
	unix('~/bin/RepeatModeler-2.0.2a/RepeatModeler -database ' + genome_short + ' -engine ncbi -pa 10', shell=True)

	#Combine databases
	os.chdir('../../CustomLib')
	os.mkdir(genome_short)
	os.chdir(genome_short)

	unix('cp ../../RepeatModeler/' + genome_short + '/RM_*/consensi.fa.classified RepeatModeler_' + genome_short + '_consensi.fa.classified.lib', shell=True)
	unix('cp ../../data/RepeatMaskerInsecta.lib .', shell=True)
	unix('cp ../../data/SINEs.bnk .', shell=True)
	unix('cat * > CombinedLibrary.lib', shell=True)


	#Filter combined library and run RepeatMasker
	unix('seqtk seq -L 50 CombinedLibrary.lib > CombinedLibrary.lib.minlen50', shell=True)
	unix('~/bin/usearch -cluster_fast CombinedLibrary.lib.minlen50 -id 0.8 -consout CombinedLibrary.lib.minlen50.nr', shell=True)

	unix('~/bin/RepeatModeler-2.0.2a/RepeatClassifier -consensi CombinedLibrary.lib.minlen50.nr', shell=True)


	#Filter for potential protein coding genes
	os.chdir('../../protein_coding_filter')
	os.mkdir(genome_short)
	os.chdir(genome_short)

	unix('cp ../../CustomLib/' + genome_short + '/CombinedLibrary.lib.minlen50.nr.classified .', shell=True)
	unix('grep "Unknown" CombinedLibrary.lib.minlen50.nr.classified > UnknownIdslist.txt', shell=True)
	unix("sed -i 's/>//g' UnknownIdslist.txt", shell=True)	
	unix('python ~/TE_annotation/data/faSomeRecords/faSomeRecords.py --fasta CombinedLibrary.lib.minlen50.nr.classified --list UnknownIdslist.txt --outfile UnknownRepeats.fasta', shell=True)
	unix('blastx -query UnknownRepeats.fasta -db ~/TE_annotation/data/UniprotSeqs/uniprot-reviewed_yes+taxonomy_50557.fasta -evalue 1e-10 -num_threads 4 -max_target_seqs 1 -outfmt "6 qseqid sseqid evalue bitscore sgi sacc stitle" -out Blast_out.txt', shell=True)
	unix("awk -F '\t' '{print $1,$7}' Blast_out.txt  | sort | uniq | grep -i -v 'transposon' | grep -i -v 'Copia protein' | grep -i -v 'mobile element' | grep -i -v 'transposable'  | grep -i -v 'transposase' | awk '{print $1}' > Unknowns_with_Port_hit.txt", shell=True)
	unix('python ~/TE_annotation/data/faSomeRecords/faSomeRecords.py --fasta CombinedLibrary.lib.minlen50.nr.classified --list Unknowns_with_Port_hit.txt --outfile ' + genome_short + '_Rep_CombinedLibrary.lib.minlen50.nr.classified.filtered.fa --exclude', shell=True)

	#Run RepeatMasker to identify Repeats
	os.mkdir('../../repeatMasker/' + genome_short)
	unix('~/bin/RepeatMasker/RepeatMasker -pa 10 -lib ' + genome_short + '_Rep_CombinedLibrary.lib.minlen50.nr.classified.filtered.fa ~/TE_annotation/data/genomes/' + genome + ' -gff -dir ../../repeatMasker/' + genome_short + '/', shell=True)

sp_list = []
for species in glob.glob('~/TE_annotation/data/genomes/*fasta'):
	species = species.split('/')[-1]
	sp_list.append(species)


Parallel(n_jobs=10)(delayed(TE_analysis)(genome) for genome in sp_list)






