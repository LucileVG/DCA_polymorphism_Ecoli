import sys
import os
import multiprocessing

from pathlib import Path
from joblib import Parallel, delayed
src_dir = os.path.join(os.getcwd(), 'src')
sys.path.append(src_dir)
from config import (cl_dv_dca_msa_path, cl_dv_dca_models_path, tmp_path,
cl_dv_tmp_analysis_folder, cl_dv_homologs_path, results_folder)
from utils import (process_gene, full_seq, merge_full_seq_dfs,
 couplings, merge_couplings_dfs)

def main():
    num_cores = multiprocessing.cpu_count()
    Path(results_folder).mkdir(parents=True, exist_ok=True)
    msa_folder = Path(cl_dv_dca_msa_path)
    aln_folder = Path(cl_dv_homologs_path)
    dca_path = Path(cl_dv_dca_models_path)
    output_folder = Path(tmp_path) / cl_dv_tmp_analysis_folder
    output_folder.mkdir(parents=True, exist_ok=True)
    genes = [str(file.stem).split(".")[0] for file in dca_path.iterdir() if (msa_folder / f"{file.stem}").is_file()]
    Parallel(n_jobs=num_cores)(delayed(process_gene)(gene, output_folder, msa_folder, aln_folder, dca_path) for gene in genes)
    folders = [folder for folder in output_folder.iterdir() if (folder / "sequence_scores.csv").is_file()]
    Parallel(n_jobs=num_cores)(delayed(full_seq)(gene) for gene in folders)
    merge_full_seq_dfs(output_folder, Path(f"{results_folder}/full_seq_single_muts.csv"))
    Parallel(n_jobs=num_cores)(delayed(couplings)(gene, dca_path, aln_folder) for gene in folders)
    merge_couplings_dfs(output_folder, Path(f"{results_folder}/couplings.csv"))

if __name__ == "__main__":
    main() 