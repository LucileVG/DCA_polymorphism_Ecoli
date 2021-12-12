# Global variables
valid_aa = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','-']
AA_list = valid_aa + ['Other']

Codons_list = list()
for i in list("ATGC"):
    for j in list("ATGC"):
        for k in list("ATGC"):
            Codons_list.append(i+j+k)
Codons_list += ["---","Other"]

TRANSLATION_TABLE = { 
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'}

d_aa_num= {a:i for i,a in enumerate(valid_aa)}

# Paths to folders and files

reference = "ESC_GA4805AA"

# Input data
path_reference = f"./datasets/{reference}.fasta"
local_dca_msa_path = './datasets/DCA_training_MSAs/local_strains'
cl_dv_dca_msa_path = './datasets/DCA_training_MSAs/closely_related_species'
local_dca_models_path = './datasets/DCA_models/local_strains'
cl_dv_dca_models_path = './datasets/DCA_models/closely_related_species'
local_homologs_path = f'./datasets/homologous_sequences/local_strains'
cl_dv_homologs_path = f'./datasets/homologous_sequences/closely_related_species'


# Output files

# Directories for temporary files
tmp_path = './tmp'
local_tmp_analysis_folder = "local_strains" 
cl_dv_tmp_analysis_folder = "closely_related_species"

# Directory for results
results_folder = "./results"
# stats file
stats_file_path = f'{results_folder}/stats_{reference}.csv'
# simulations file
simulations_file_path = f"{results_folder}/simulated_sites_{reference}.csv"
# mutants sites file
mutants_file_path = f"{results_folder}/mutants_sites_{reference}.csv"


# Path to program

path_julia = "julia"


