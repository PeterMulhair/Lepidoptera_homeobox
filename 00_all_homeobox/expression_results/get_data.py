import os
from subprocess import call as unix
from joblib import Parallel, delayed

os.mkdir('raw')

def download_fastq(fastq):
    os.chdir('raw')
    unix('sudo wget ' + fastq, shell = True)
    os.chdir('../')
    
ID_sp = {}
fastq_list = []        
with open('species_data.tsv') as f:
    for line in f:
        lines = line.split('\t')
        sp = lines[0]
        sp = sp.replace(' ', '_')
        data = lines[-1].strip()
        data2 = data.split('_1.fastq')[0] + '_2.fastq.gz'
        ERRID = lines[2]
        fastq_list.append(data)
        fastq_list.append(data2)
        ID_sp[ERRID] = sp

Parallel(n_jobs=20)(delayed(download_fastq)(download) for download in fastq_list)
        
for k, v in ID_sp.items():
    unix('mv ' + k + '_1.fastq.gz ' + v + '_1.fastq.gz', shell=True)
    unix('mv ' + k + '_2.fastq.gz ' + v + '_2.fastq.gz', shell=True)
