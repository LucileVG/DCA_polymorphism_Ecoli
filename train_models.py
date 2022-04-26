import sys
import os
import multiprocessing
from pathlib import Path
from joblib import Parallel, delayed

src_dir = os.path.join(os.getcwd(), 'src')
sys.path.append(src_dir)
from config import (local_dca_msa_path, cl_dv_dca_msa_path,
local_dca_models_path, cl_dv_dca_models_path, path_julia)

def run_ind():
    os.system(f"{path_julia} ./compute_freq_julia.jl")

def run_dca(path_msa, path_dca, num_cores):
    os.system(f"{path_julia} ./plmdca.jl {path_msa} {path_dca} --threads {num_cores}")

def main():
    num_cores = multiprocessing.cpu_count()
    local_dca_models_folder = Path(local_dca_models_path)
    cl_dv_dca_models_folder = Path(cl_dv_dca_models_path)
    local_dca_models_folder.mkdir(parents=True, exist_ok=True)
    cl_dv_dca_models_folder.mkdir(parents=True, exist_ok=True)
    local_dca_msa_folder = Path(local_dca_msa_path)
    cl_dv_dca_msa_folder = Path(cl_dv_dca_msa_path)
    run_ind()
    for input_msa in local_dca_msa_folder.iterdir():
        run_dca(input_msa,  local_dca_models_folder / f"{input_msa.name}.npz", num_cores)
    for input_msa in cl_dv_dca_msa_folder.iterdir():
        run_dca(input_msa, cl_dv_dca_models_folder / f"{input_msa.name}.npz", num_cores)

if __name__ == "__main__":
    main() 
