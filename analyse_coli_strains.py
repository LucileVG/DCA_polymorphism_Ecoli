import sys
import os
import logging
import pandas as pd
from pathlib import Path
from joblib import Parallel, delayed
src_dir = os.path.join(os.getcwd(), 'src')
sys.path.append(src_dir)
from utils import (fasta2dict, read_fasta, write_fasta, 
directory_content, accessible_mutations, accessible_codons, 
synonymous_mutations,  homologs_analysis, distant_homologs_analysis, 
gene_single_mutant_scores, gather_possible_mutants, 
compute_site_stats, merge_dfs, gather_info_sites, double_muts, 
merge_df,IPR)
from config import (name_ref_proteome, path_ref_proteome, tmp_path, 
msa_path_hhblits, path_dca_par, dna_msa_path, msa_tmp_analysis_folder, 
stats_file_path, simulations_file_path, mutants_file_path, make_folder)

import multiprocessing

# Workflow

def make_gene_folders(path):
    '''
    Create a temporary folder with subfolders at given path.
    '''
    Path(path / "homologs_analysis" / "sequences").mkdir(parents=True, exist_ok=True)
    Path(path / "reference_analysis").mkdir(parents=True, exist_ok=True)
    Path(path / "distant_homologs").mkdir(parents=True, exist_ok=True)
    Path(path / "stats").mkdir(parents=True, exist_ok=True)


def test_dca_model(path_to_file):
    if(path_to_file.is_file()):
        try:
            matrices = np.load(path_to_file)
            return True
        except:
            return False
    return False


def main():
    num_cores = 15 #multiprocessing.cpu_count()
    threshold = 0.05
    make_folder()
    homologs = directory_content(Path(dna_msa_path))
    reference_file = "dna_{}".format(name_ref_proteome)
    reference = fasta2dict(Path(path_ref_proteome) / reference_file)
    logging.info('Initialize temporary folders for MSAs analysis')
    #genes = list()
    #for file in homologs:
    #    gene = file.stem
    #    dca_model = Path(os.path.join(path_dca_par, gene+".fa.npz"))
    #    distant_homolog_file = "{}.fa".format(gene)
    #    if(gene in reference.keys() and test_dca_model(dca_model) and (Path(msa_path_hhblits) / distant_homolog_file).is_file()):
    #        gene_folder =  Path(tmp_path) / msa_tmp_analysis_folder / gene
    #        genes.append(gene_folder)
    #        make_gene_folders(gene_folder)
    #        write_fasta(gene_folder / "reference_analysis" / "reference_sequence.fasta", [reference[gene]])
    #        write_fasta(gene_folder / "homologs_analysis" / "homologs.fasta", read_fasta(file))
    #        write_fasta(gene_folder / "distant_homologs" / "distant_homologs.fasta", read_fasta(Path(msa_path_hhblits) / distant_homolog_file))
    
    genes = [folder for folder in (Path(tmp_path) / msa_tmp_analysis_folder).iterdir()]
    
    #logging.info('Compute DCA scores')
    #Parallel(n_jobs=num_cores)(delayed(gene_single_mutant_scores)(gene) for gene in genes)
    #logging.info('Compute 1-SNP AA mutations')
    #Parallel(n_jobs=num_cores)(delayed(accessible_mutations)(gene / "reference_analysis") for gene in genes)
    #logging.info('Compute 1-SNP codons')
    #Parallel(n_jobs=num_cores)(delayed(accessible_codons)(gene / "reference_analysis") for gene in genes)
    #logging.info('Compute synonymous mutations')
    #Parallel(n_jobs=num_cores)(delayed(synonymous_mutations)(gene / "reference_analysis") for gene in genes)
    #logging.info('Alignment and analysis of homologs')
    #Parallel(n_jobs=num_cores)(delayed(homologs_analysis)(gene) for gene in genes)
    #logging.info('Analysis of distant homologs')
    #Parallel(n_jobs=num_cores)(delayed(distant_homologs_analysis)(gene / "distant_homologs") for gene in genes)
    #logging.info('Stats computations')
    #Parallel(n_jobs=num_cores)(delayed(gather_possible_mutants)(gene) for gene in genes)
    #dfs = [pd.read_csv(gene / "stats" / "mutants_sites.csv") for gene in genes]
    #pd.concat(dfs).to_csv(mutants_file_path, index=None)
    #Parallel(n_jobs=num_cores)(delayed(compute_site_stats)(gene, threshold) for gene in genes)
    #dfs = [merge_dfs(gene, threshold) for gene in genes]
    #pd.concat(dfs).to_csv(stats_file_path, index=None)
    #logging.info('Parallel gathering of information for simulations')
    #Parallel(n_jobs=num_cores)(delayed(gather_info_sites)(gene_folder) for gene_folder in genes)
    #dfs = [pd.read_csv(gene / "stats" / "simulated_sites.csv") for gene in genes]
    #pd.concat(dfs).to_csv(simulations_file_path, index=None)
    
    input_folder = Path("./tmp/msa_analysis")
    output_folder = Path("./data/local")
    folders = [folder for folder in input_folder.iterdir() if (folder / "reference_analysis" / "DCA_proba.csv").is_file()]
    #Parallel(n_jobs=num_cores)(delayed(double_muts)(gene) for gene in folders)
    merge_df("double_mut_epistasis.csv", input_folder, output_folder)
    
    Parallel(n_jobs=num_cores)(delayed(IPR)(gene) for gene in folders)
    merge_df("IPR.csv", input_folder, output_folder)


if __name__ == "__main__":
    main() 