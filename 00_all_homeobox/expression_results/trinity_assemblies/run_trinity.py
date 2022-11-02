import glob
from subprocess import call as unix
from joblib import Parallel, delayed

#sp_list = ['Biston_betularia', 'Limenitis_camilla', 'Nymphalis_urticae', 'Pararge_aegeria']#, ['Pieris_rapae', 'Vanessa_atalanta', 'Vanessa_cardui']
sp_list = ['Pieris_rapae', 'Vanessa_atalanta', 'Vanessa_cardui']  

def trinity(species):
    print(species)
    unix('Trinity --seqType fq --max_memory 50G --left ../trim_map_reads/' + species + '_1_paired_trim_reform.fastq --right ../trim_map_reads/' + species + '_2_paired_trim_reform.fastq --CPU 2 --output /home/zoo/zool2500/tree_of_life/trinity_out_dir_' + species, shell=True)
    print(species, 'complete')
    
Parallel(n_jobs=4)(delayed(trinity)(sp) for sp in sp_list)
