import sys
import os
import multiprocessing

from pathlib import Path
from joblib import Parallel, delayed
src_dir = os.path.join(os.getcwd(), 'src')
sys.path.append(src_dir)
#from config import *
from utils import process_gene, full_seq, merge_full_seq_dfs, couplings, merge_couplings_dfs

def main():
    num_cores = 15 #multiprocessing.cpu_count()
    msa_folder = Path("./data/distant_homologs_2/msa")
    aln_folder = Path("./data/local/aln_on_ref/")
    dca_path = Path("./data/distant_homologs_2/DCAparameters")
    output_folder = Path("./tmp/Distant_species")
    output_folder.mkdir(parents=True, exist_ok=True)
    genes = [str(file.stem).split(".")[0] for file in dca_path.iterdir() if(len(str(file.stem).split(".")[0])>2 and str(file.stem)[:2]=="PF") and (msa_folder / f"{str(file.stem).split('.')[0]}.fasta").is_file()]
    Parallel(n_jobs=num_cores)(delayed(process_gene)(gene, output_folder, msa_folder, aln_folder, dca_path) for gene in genes)
    folders = [folder for folder in output_folder.iterdir() if (folder / "sequence_scores.csv").is_file()]
    Parallel(n_jobs=num_cores)(delayed(full_seq)(gene) for gene in folders)
    merge_full_seq_dfs(output_folder, Path("./data/local/full_seq_single_muts.csv"))
    Parallel(n_jobs=num_cores)(delayed(couplings)(gene, dca_path, aln_folder) for gene in folders)
    merge_couplings_dfs(output_folder, Path("./data/local/couplings.csv"))

if __name__ == "__main__":
    main() 