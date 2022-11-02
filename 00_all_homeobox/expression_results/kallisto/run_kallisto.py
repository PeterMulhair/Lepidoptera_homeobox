import glob
from subprocess import call as unix
from joblib import Parallel, delayed

#sp_list = ['Biston_betularia', 'Limenitis_camilla', 'Nymphalis_urticae', 'Pararge_aegeria']#, ['Pieris_rapae', 'Vanessa_atalanta', 'Vanessa_cardui']
sp_list = ['Pieris_rapae', 'Vanessa_atalanta', 'Vanessa_cardui']  

def kallisto(species):
    print(species)
    #unix('kallisto index ../trinity_assemblies/trinity_out_dir_' + species + '/Trinity.fasta -i /home/zoo/zool2500/tree_of_life/' + species + '_transcriptome', shell=True)
    unix('kallisto quant -i ' + species + '_transcriptome -o /home/zoo/zool2500/tree_of_life/' + species + '_output --threads 10 ../trim_map_reads/' + species + '_1_paired_trim_reform.fastq ../trim_map_reads/' + species + '_2_paired_trim_reform.fastq', shell=True)
       
    
Parallel(n_jobs=4)(delayed(kallisto)(sp) for sp in sp_list)
