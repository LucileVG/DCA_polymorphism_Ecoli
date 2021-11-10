import numpy as np
from pathlib import Path
import os

#aa list and dictionaries
valid_aa = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','-']
aa_3= ['ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR','-']

AA_list = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','-','Other']

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
d_num_aa= {i:a for i,a in enumerate(valid_aa)}
d_3to1 = {a3:a1 for a3,a1 in zip(aa_3,valid_aa)}

def get_list_ref( list_genes = 'list_genes_more61000.txt'):
    list_ref_prot = []
    #list ref protein to consider! (ONLY CORE GENOME!)
    f = open(os.path.join(path_ref_proteome,list_genes), "r")
    for line in f.readlines():
        list_ref_prot.append("ref_"+line.rstrip())
    f.close()
    return list_ref_prot


#SET ALL PATHs
name_ref_proteome = "ESC_GA4805AA.fa"
#name_ref_proteome = "ESC_consensus.fa"

path_ref_proteome = "./data/refstrain_"+name_ref_proteome[:-3]
tmp_path = './tmp'  #make directory for tmp files (reference sequences, a3m file from hhblits)
#for hhblits
path_db_hhblits = "/media/data/db/UniRef30_2020_03/UniRef30_2020_03"
name_db_hhblits = "UniRef30_2020_03"
path_hhblits = "/home/admin/Documents/programs/hh-suite/build/bin/hhblits"
msa_path_hhblits = './data/distant_homologs/msa'
#for dca
path_dca_par = './data/distant_homologs/DCAparameters'
#for phmmer
path_db_phmmer= "/media/data/db/StrainEcoRef100/all_strains_proteomes_uniq_with_cluster_members.fa"
path_phmmer = "/home/admin/Documents/programs/hmmer-3.3.1/src/phmmer"
#msa for energy-computation (phmmer or random)
msa_path_phmmer = './data/local/msa_'+name_ref_proteome[:-3]
msa_path_ind = './data/local/msa_'+name_ref_proteome[:-3]+'_ind'
msa_path_random = './data/local/msa_'+name_ref_proteome[:-3]+'_random'
#folder containing enegy data
e_path = './data/e_'+name_ref_proteome[:-3]
#folder containing entropy data
s_path = './data/s_'+name_ref_proteome[:-3]

#folder with dna msa of local homologs
dna_msa_path = './data/local/dna_msa_'+name_ref_proteome[:-3]
#temporary folder for gene analysis
msa_tmp_analysis_folder = "msa_analysis" 
#stats file
stats_file_path = './data/local/stats_{}.csv'.format(name_ref_proteome[:-3])
#simulations file
simulations_file_path = "./data/local/simulated_sites_{}.csv".format(name_ref_proteome[:-3])
#mutants sites file
mutants_file_path ="./data/local/mutants_sites_{}.csv".format(name_ref_proteome[:-3])
#supplementary figures folder
supp_figures_path = "./imgs/supplementary_figures"

def make_folder():
    #distant
    Path(tmp_path).mkdir(parents=True, exist_ok=True)
    Path(msa_path_hhblits).mkdir(parents=True, exist_ok=True)
    Path(path_dca_par).mkdir(parents=True, exist_ok=True)
    #local
    Path(msa_path_phmmer).mkdir(parents=True, exist_ok=True)
    Path(msa_path_ind).mkdir(parents=True, exist_ok=True)
    Path(msa_path_random).mkdir(parents=True, exist_ok=True)
    Path(e_path).mkdir(parents=True, exist_ok=True)
    Path(s_path).mkdir(parents=True, exist_ok=True)
    return 0


