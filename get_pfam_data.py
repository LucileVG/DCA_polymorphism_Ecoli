import pickle
import numpy as np
#from prody import *
#from matplotlib.pylab import *
from Bio import AlignIO
import re
from Bio import SeqIO
import sys
from collections import OrderedDict
import os
import pandas as pd
#import pyhmmer
#from hmm_profile import reader
#from hmm_profile import writer
import subprocess
from joblib import Parallel, delayed
import multiprocessing
import os
import time


####################################################################################################
# PATH TO EXECUTABLES
#set path to hmmalign
path_hmmalign = '/work/FAC/FBM/LLB/dgfeller/epitope_pred/gcroce1/dca/hmmer-3.3.2/src/hmmalign'
#set path to julia bin file
path_julia = 'julia'
#set path to reference E. coli strain (main reference) -> avoid overfitting
path_coli = './ESC_GA4805AA.fa'
#set path to folder containing reference MSA (Ecoli + distant species) -> avoid overfitting also for distant species
path_ref_folder = "./Aln_pfam"
#set folder for temporary files
tmp_folder = './tmp'
#set filtering (min sequences dist from reference sequences) -> avoid overfitting
hole_perc = 0.10
#set num cores (for parallele computing -> automatically select all the CPUs)
num_cores = multiprocessing.cpu_count()

####################################################################################################
# FUNCTIONS
remove_lower = lambda text: re.sub('[a-z]', '', text)

def compute_dist(ref_seq, seq):
    distance = sum([1 for x, y in zip(ref_seq, seq) if x.lower() != y.lower()])
    return distance

def compute_seqid(ref_seq, seq):
    '''return the sequence identity (seqid) '''
    distance = compute_dist(ref_seq,seq)
    distance /= len(seq)
    seqid = 1 - distance
    return seqid

def adapt_pfam(path_pfam_stk, seq_aligned, list_ref_id, path_pfam_out, mask):
    ''' read pfam stk, remove seq in seq_aligned and close-by seq (10%)
        mask: to rm positions that are gapped in the reference '''
    removed_seq = 0
    #mask the reference sequence (and store masked ref seq)
    prot_name = path_pfam_stk.split("_")[0] + "_" + path_pfam_stk.split("_")[1]
    pfam_name = path_pfam_stk.split("_")[2]
    ref_masked_coli = open("{0}-{1}.fasta".format(prot_name, pfam_name), "w")
    seq_aligned_masked = []
    print(list_ref_id)
    for ref_id, ref_seq in zip(list_ref_id, seq_aligned):
        masked_ref_seq  = ''.join(ref_seq[x] for x in range(len(ref_seq)) if x not in mask)
        seq_aligned_masked.append(masked_ref_seq)
        print(">{0}".format(ref_id), file = ref_masked_coli)
        print(masked_ref_seq, file = ref_masked_coli)
    ref_masked_coli.close()
    #now adapt pfam
    f_out = open(path_pfam_out, 'w')
    with open(path_pfam_stk) as fp:
        lines = fp.readlines()
        for line in lines:
            #skip  comments
            if line[0] == "#" or line[0] == "/":
                continue
            id_seq = line.split()[0]
            seq = line.split()[-1]
            #rm '.' and lower chars
            seq = remove_lower(seq.replace(".", ""))
            #here mask sequence
            masked_seq = ''.join(seq[x] for x in range(len(seq)) if x not in mask)
            #filter -> rm seq close to the ref_seq (masked)
            ok_seq = True
            for ref_seq in seq_aligned_masked:
                seqid = compute_seqid(masked_ref_seq, masked_seq)
                if seqid < hole_perc:
                    ok_seq = False
            if ok_seq:
                f_out.write('>{0}\n{1}\n'.format(id_seq, masked_seq))
            else:
                removed_seq += 1
    f_out.close()
    print('removed {0} seq'.format(removed_seq))
    return 0

def get_hmm_file(prot, pfam, path_pfam_hmm):
    ''' download hmm pfam '''
    if os.path.exists('./tmp/hmm_{0}_{1}'.format(prot, pfam)):
        return 0
    #to avoid overlap when parallel
    subprocess.run(['curl', '-o', './tmp/hmm_{0}_{1}'.format(prot, pfam), path_pfam_hmm])
    return 0

def align_ref_pfam(list_ref_seq, path_hmm, pfam):
    '''align reference sequences to pfam hmm'''
    list_seq_aligned = []
    for seq in list_ref_seq:
        #write tmp seq
        path_tmp_ref_seq = os.path.join(tmp_folder, 'ref_seq.fa')
        print(">{0}\n{1}\n".format("seq", seq), file = open(path_tmp_ref_seq, 'w'))
        # run hmmalign
        path_out_hmm = os.path.join(tmp_folder, "tmp_hmmalign_{0}.txt".format(pfam))
        FOUT = open(path_out_hmm, 'w')
        subprocess.run([path_hmmalign, path_hmm, path_tmp_ref_seq], stdout=FOUT, stderr=subprocess.STDOUT)
        FOUT.close()
        # parse hmm file
        alignment = AlignIO.read(open(path_out_hmm), "stockholm")
        for record in alignment:
            #print(record.id)
            seq_aligned_to_msa  = remove_lower(str(record.seq))
            #print(seq_aligned_to_msa, len(seq_aligned_to_msa))
            list_seq_aligned.append(seq_aligned_to_msa)
    return list_seq_aligned


#### DCA ANALYSIS ####
def run_dca(path_msa, path_dca, num_cores):
    subprocess.run([path_julia, '--threads', str(num_cores), './plmdca.jl', path_msa, path_dca])
    #'''run plmDCA, python wrapper to plmDCA_julia (Pagnani version )'''
    ##julia.install(path_julia)
    #from julia import Julia
    #from julia.api import Julia
    ##compiled_modules == True doesn't work.. I need to recompile the module each time (slow)
    #jl = Julia(runtime=path_julia, compiled_modules=False)
    #from julia import PlmDCA
    #dca_res = PlmDCA.plmdca(path_msa, theta = 0.2, verbose = False)
    #J = dca_res.Jtensor
    #h = dca_res.htensor
    #fn_dca = dca_res.score
    #return h,J, fn_dca
    return 0

def load_dca(prot, pfam):
    ''' load the dca model in ./pfam_dca, return h,J  '''
    loaded = np.load('./pfam_dca/{0}-{1}.fa.npz'.format(prot, pfam))
    h, J =loaded['h_dca'], loaded['J_dca']
    return h, J

def get_mask(seq_aligned):
    ''' get position that are gapped in the reference '''
    idx_gap = []
    for seq in seq_aligned:
        #print(seq)
        for idx, char in enumerate(list(seq)):
            if char == "-":
                idx_gap.append(idx)
        #print(idx_gap)
    #check all reference have gap in same positions!
    mask, idx_gap_counts = np.unique(idx_gap, return_counts = True)
    #if len(np.unique(idx_gap_counts)) > 0:
    #    assert len(np.unique(idx_gap_counts)) == 1
    return mask


#idea: read coli file (ref_coli), get also ref_yers
def compute_pfam_dca(keys, only_ecoli_as_ref = False):
    prot, pfam = keys
    #set ref seq (
    #ref_coli = d_coli[prot]
    #list_ref_seq = [str(ref_coli)]#, str(ref_yers)]
    #add lucile files (for multihole)
    list_ref_seq = []
    list_ref_id = []
    name_file = "{0}-{1}.fasta".format(pfam,prot)
    path_ref = os.path.join(path_ref_folder, name_file)
    if only_ecoli_as_ref == True:
        # this is to keep only ECOLI ESC_GA4805AA  as reference genome
        with open(path_ref) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                list_ref_id.append(record.id)
                seq = record.seq
                list_ref_seq.append(seq)
                break # this is to keep only ECOLI ESC_GA4805AA
    if only_ecoli_as_ref == False:
        # this is to ECOLI ESC_GA4805AA + distant species as reference genome
        with open(path_ref) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                list_ref_id.append(record.id)
                seq = record.seq
                list_ref_seq.append(seq)
    print('0. Download pfam: {0}'.format(pfam))
    #download
    path_pfam_sth= 'http://pfam.xfam.org/family/{0}/alignment/full/gzipped'.format(pfam)
    subprocess.run(['curl', '-o', pfam+'_full.sth.gz', path_pfam_sth ])
    ##gunzip
    subprocess.run(['gunzip', pfam+'_full.sth.gz'])
    path_pfam_stk= '{0}_{1}_full.sth'.format(prot, pfam)
    ##final name
    subprocess.run(['mv', pfam+"_full.sth", path_pfam_stk])
    print('1. Align refs to pfam hmm')
    path_pfam_hmm = 'http://pfam.xfam.org/family/{0}/hmm'.format(pfam)
    get_hmm_file(prot, pfam, path_pfam_hmm)
    seq_aligned = align_ref_pfam(list_ref_seq, './tmp/hmm_{0}_{1}'.format(prot, pfam), pfam)
    mask = get_mask(seq_aligned)
    print('2. filter pfam msa')
    path_pfam_filter = './pfam_msa/{0}-{1}.fa'.format(prot, pfam)
    adapt_pfam(path_pfam_stk, seq_aligned, list_ref_id, path_pfam_filter, mask)
    #rm useless files
    subprocess.run(['rm',path_pfam_stk])
    print('3. Train DCA model')
    path_julia_cluster = "/dcsrsoft/spack/hetre/v1.1/spack/opt/spack/linux-rhel8-zen2/gcc-9.3.0/julia-1.6.0-jjb27b3skndm6fwbecdzgdxpd7piudje/bin/julia"
    path_out_dca= './pfam_dca/{0}-{1}.fa.npz'.format(prot, pfam)
    run_dca(path_pfam_filter, path_out_dca,num_cores)
    ##load dca model
    ##a, b = load_dca(prot, pfam)
    #except:
    #    print("ERROR with {0}".format(keys))
    return 0

if __name__== "__main__":
    prot_input = sys.argv[1]
    pfam_input = sys.argv[2]
    print(prot_input, pfam_input)
    #0. read reference (coli, yers)
    #make dict d[(prot, pfam)] = seq
    d_coli = {}
    with open(path_coli) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seq = record.seq
            if record.id.split(" ")[0][4:] == prot_input:
                d_coli[prot_input] = seq
    #main function -> download pfam -> filter -> train dca model
    compute_pfam_dca((prot_input, pfam_input), only_ecoli_as_ref = False)
