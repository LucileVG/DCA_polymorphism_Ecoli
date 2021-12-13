import os
import numpy as np
import pandas as pd

from Bio import SeqIO
from pathlib import Path
from collections import defaultdict
from config import valid_aa, d_aa_num, Codons_list, TRANSLATION_TABLE, AA_list

#from config import *

#import sys
#import random
#import seaborn as sns
#import matplotlib.pyplot as plt
#import gzip
#import subprocess
#import tarfile
#from io import BytesIO
#from joblib import Parallel, delayed
#import multiprocessing
#from tqdm import tqdm

def directory_content(path):
    '''
    Return a list of paths to all the files contained in a directory.
    '''
    return [x for x in path.iterdir()]

def fasta2dict(path):
    '''
    Read a fasta file and returns its sequences as a dict of SeqRecord objects.
    '''
    with open(path, "r") as input:
        return {record.id:record for record in SeqIO.parse(input, "fasta")}

def write_fasta(path,sequences):
    '''
    Write a fasta file with the sequences of the list sequences in the same order.
    '''
    with open(path, "w") as output:
        SeqIO.write(sequences,output,"fasta")

def read_fasta(path):
    '''
    Read a fasta file and returns its sequences as a list of SeqRecord objects.
    '''
    with open(path, "r") as input:
        return list(SeqIO.parse(input, "fasta"))

def read_dca_par(dca_models_path, ref_prot):
    matrices = np.load(os.path.join(dca_models_path, ref_prot+".fa.npz"))
    return matrices["h_dca"], matrices["J_dca"]

def translate(sequence):
    '''
    Translate the dna sequence to amino acid sequence. Returns None if its length is not a multiple of 3.
    '''
    Translation_table = { 
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
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
    '---':'-'}
    if(len(sequence)%3==0):
        upper_sequence = sequence.upper()
        translated_sequence = ""
        for i in range(len(upper_sequence)//3):
            codon = upper_sequence[i*3:(i+1)*3]
            aa = "X"
            if(codon in Translation_table.keys()):
                aa = Translation_table[codon]
            translated_sequence = translated_sequence+aa
        return translated_sequence
    return None

def all_standard_aa(seq):
    '''Return True if sequence contains only standard amino acids.'''
    for char in seq:
        if((char not in valid_aa) and char !='-'):
            return False
    return True

def compute_energy(seq, h, J):
    '''
    Compute energy of the given sequence.
    '''
    if all_standard_aa(seq):
        E = 0
        for idx_aa1 in range(0, len(seq)):
            aa1 = seq[idx_aa1]
            E -= h[d_aa_num[aa1], idx_aa1]
            for idx_aa2 in range(idx_aa1+1, len(seq)):
                aa2 = seq[idx_aa2]
                E -= J[d_aa_num[aa1], d_aa_num[aa2], idx_aa1, idx_aa2]
        return E

def perform_DCA(sequence, J, h):
    '''
    Compute energy for any single mutant and store results in a matrix.
    '''
    n = len(sequence)
    energy_ref = compute_energy(sequence, h, J)
    energy_matrix = np.zeros((n,len(valid_aa)))
    for i in range(n):
        for j in range(len(valid_aa)):
            if(sequence[i]!=valid_aa[j]):
                mutant = sequence[:i]+valid_aa[j]+sequence[i+1:]
                energy_matrix[i,j] = compute_energy(mutant, h, J)-energy_ref
    return energy_matrix

def entropy(x):
    '''
    Compute single entropy term.
    '''
    if(x>0):
        return -np.log(x)/np.log(2)*x
    return 0

def matrix2df(matrix, sequence, output_folder):
    '''
    Compute quantities from DCA scores and write results in corresponding csv files.
    '''
    df = pd.DataFrame(matrix)
    df.columns = valid_aa
    df["Context"] = list(sequence)
    for character in valid_aa[:-1]:
        df["{}_exp".format(character)] = np.exp(-df[character])
    exp_columns = ["{}_exp".format(character) for character in valid_aa[:-1]]
    df["Sum"] = df[exp_columns].apply(sum,axis=1)
    for character in valid_aa[:-1]:
        df["{}_proba".format(character)] = df["{}_exp".format(character)]/df["Sum"]
        df["{}_entropy".format(character)] = df["{}_proba".format(character)].apply(entropy)
    proba_columns = ["{}_proba".format(character) for character in valid_aa[:-1]]
    entropy_columns = ["{}_entropy".format(character) for character in valid_aa[:-1]]
    df["CD_Entropy"] = df[entropy_columns].apply(sum,axis=1) 
    df[["Context"]+valid_aa].to_csv(output_folder / "DCA_scores.csv",index=None)
    df[["Context"]+proba_columns].to_csv(output_folder / "DCA_proba.csv",index=None)
    df[["Context","CD_Entropy"]].to_csv(output_folder / "CD_Entropy.csv",index=None)

def gene_single_mutant_scores(dca_models_path, output_folder):
    '''
    For a given gene, compute DCA scores and derivated quantities and write results in corresponding csv files.
    '''
    h, J = read_dca_par(dca_models_path, output_folder.name)
    output_directory = output_folder / "reference_analysis"
    sequence = read_fasta(output_directory / "reference_sequence.fasta")[0]
    AA_sequence = translate(str(sequence.seq))
    energy_matrix = perform_DCA(AA_sequence, J, h)
    matrix2df(energy_matrix, AA_sequence, output_directory)

def accessible_mutations(output_folder):
    '''
    Generate a dataframe with accessible AA from reference sequence.
    '''
    columns = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','*']
    sequence = str(read_fasta(output_folder / "reference_sequence.fasta")[0].seq)
    upper_sequence = sequence.upper()
    n = len(upper_sequence)//3
    results = pd.DataFrame(np.zeros((n,21)))
    results.columns = columns
    results["Codon"] = ""
    for i in range(n):
        codon = upper_sequence[i*3:(i+1)*3]
        results.at[i,"Codon"] = codon
        if(codon in TRANSLATION_TABLE.keys()):
            mutated_codons = list()
            for character in list("ACGT"):
                mutated_codons.append(character + codon[1:])
                mutated_codons.append(codon[0] + character + codon[2])
                mutated_codons.append(codon[:2] + character)
            for new_codon in mutated_codons:
                results.loc[i,TRANSLATION_TABLE[new_codon]] += 1
    results[["Codon"]+columns].to_csv(output_folder / "accessible_mutants.csv", index=None)

def accessible_codons(output_folder):
    '''
    Generate a dataframe with accessible codons from reference sequence.
    '''
    Codons_list = list(TRANSLATION_TABLE.keys())
    sequence = str(read_fasta(output_folder / "reference_sequence.fasta")[0].seq)
    upper_sequence = sequence.upper()
    n = len(upper_sequence)//3
    m = len(Codons_list)
    results = pd.DataFrame(np.zeros((n,m)))
    results.columns = Codons_list
    results["Codon"] = ""
    for i in range(n):
        codon = upper_sequence[i*3:(i+1)*3]
        results.at[i,"Codon"] = codon
        if(codon in TRANSLATION_TABLE.keys()):
            mutated_codons = list()
            for character in list("ACGT"):
                mutated_codons.append(character + codon[1:])
                mutated_codons.append(codon[0] + character + codon[2])
                mutated_codons.append(codon[:2] + character)
            for new_codon in mutated_codons:
                results.at[i,new_codon] = 1
    results[["Codon"]+Codons_list].to_csv(output_folder / "accessible_codons.csv", index=None)

def synonymous_mutations(output_folder):
    '''
    Generate a dataframe with synonymous mutations of reference sequence.
    '''
    Codons_list = list(TRANSLATION_TABLE.keys())
    sequence = str(read_fasta(output_folder / "reference_sequence.fasta")[0].seq)
    upper_sequence = sequence.upper()
    n = len(upper_sequence)//3
    m = len(Codons_list)
    results = pd.DataFrame(np.zeros((n,m)))
    results.columns = Codons_list
    results["Codon"] = ""
    for i in range(n):
        codon = upper_sequence[i*3:(i+1)*3]
        results.at[i,"Codon"] = codon
        if(codon in TRANSLATION_TABLE.keys()):
            aa = TRANSLATION_TABLE[codon]
            for new_codon in Codons_list:
                if(TRANSLATION_TABLE[new_codon]==aa):
                    results.at[i-1,new_codon] = 1
    results[["Codon"]+Codons_list].to_csv(output_folder / "synonymous_mutants.csv", index=None)

def extract_unique_sequences(output_folder, reference):
    '''
    Extract unique DNA sequences and write them in a fasta file with unique ids.
    A csv containing the correspondance between ids and sequence names is also created.
    '''
    unique_records = defaultdict(list)
    with open(output_folder / "homologs.fasta", "r") as input_fasta:
        for record in SeqIO.parse(input_fasta, "fasta"):
            unique_records[str(record.seq)].append(record.id)
    unique_id = 0
    refseq = str(read_fasta(reference)[0].seq)
    with open(output_folder / "sequences" / "sequence_unique_id.csv","w") as csv:
        csv.write("Name;Id\n")
        with open(output_folder / "sequences" / "unique_DNA.fasta","w") as fasta:
            fasta.write(">{}\n".format("REFERENCE"))
            fasta.write("{}\n".format(refseq.lower()))
            for key in unique_records.keys():
                for sequence_id in unique_records[key]:
                    csv.write("{};{}\n".format(sequence_id,unique_id))
                fasta.write(">{}\n".format(unique_id))
                fasta.write("{}\n".format(key))
                unique_id += 1

def translate_unique_sequences(output_folder):
    '''
    Translate a DNA fasta file to an AA fasta file.
    '''
    DNA_sequences = read_fasta(output_folder / "unique_DNA.fasta")
    AA_sequences = dict()
    for dna_seq in DNA_sequences:
        AA_sequences[dna_seq.id] = translate(str(dna_seq.seq))
    with open(output_folder / "unique_AA.fasta","w") as fasta:
        for key in AA_sequences:
            fasta.write(">{}\n".format(key))
            fasta.write("{}\n".format(AA_sequences[key]))

def align_AA_sequences(output_folder):
    '''
    Read fasta file of unique AA sequences, align them and write them in another fasta file.
    '''
    AA_hits = output_folder / "unique_AA.fasta"
    aln_AA = output_folder / "unique_AA_aln.fasta"
    os.system('mafft --quiet --retree 1 --thread 1 {} > {}'.format(AA_hits, aln_AA))

def reverse_translate(nucl_seq,aln_prot_seq):
    '''
    Reverse-align a dna sequence from an aligned protein sequence.
    '''
    start_gap = len(aln_prot_seq)-len(aln_prot_seq.lstrip('-')) # Number of gaps at the beginning of the aligned protein sequence
    end_gap = len(aln_prot_seq)-len(aln_prot_seq.rstrip('-')) # Number of gaps at the end of the aligned protein sequence
    aln_prot_seq = aln_prot_seq.lstrip('-').rstrip('-') # We remove starting and ending gaps
    values = [len(AA)*3 for AA in aln_prot_seq.split("-")] # Each element of the list is the number of aminoacids between two consecutive gaps (can be 0 if there are several gaps one after the other) multiplied by 3 in order to have the number of nucleotides
    aln_nucl_seq = "---"*start_gap # Start the aligned nucleotidic sequence
    count = 0
    for bps in values: # Loop to add nucleotides and gaps
        aln_nucl_seq = aln_nucl_seq+nucl_seq[count:count+bps]+"---"
        count += bps
    return aln_nucl_seq.rstrip('---')+"---"*end_gap # Add the ending gaps

def align_DNA_sequences(output_folder):
    '''
    Reverse align DNA sequences in a fasta file from their AA alignment.
    '''
    aln_AA_seqs = read_fasta(output_folder / "unique_AA_aln.fasta")
    raw_DNA_seqs = fasta2dict(output_folder / "unique_DNA.fasta")
    with open(output_folder / "unique_DNA_aln.fasta","w") as aln_DNA_seqs:
        for AA_seq in aln_AA_seqs:
            DNA_seq = raw_DNA_seqs[AA_seq.id]
            aln_DNA_seqs.write(">{}\n".format(AA_seq.id))
            aln_DNA_seqs.write("{}\n".format(reverse_translate(str(DNA_seq.seq),str(AA_seq.seq))))

def remove_reference(fastafile):
    '''
    Remove reference from fasta file and remove inserts in other sequences.
    '''
    aln_seq = fasta2dict(fastafile)
    reference = str(aln_seq["REFERENCE"].seq)
    indexes = [i for i in range(len(reference)) if reference[i]!='-']
    with open(fastafile,"w") as fasta:
        for key in aln_seq.keys():
            if(key != "REFERENCE"):
                fasta.write(">{}\n".format(key))
                sequence = str(aln_seq[key].seq)
                fasta.write("{}\n".format(''.join([sequence[i] for i in indexes])))

def retrieve_all_sequences(output_folder):
    '''
    Write all selected aligned DNA and AA sequences to fasta files.
    '''
    df = pd.read_csv(output_folder / "sequence_unique_id.csv",sep=";")
    DNA_sequences = fasta2dict(output_folder / "unique_DNA_aln.fasta")
    AA_sequences = fasta2dict(output_folder / "unique_AA_aln.fasta")
    with open(output_folder / "all_DNA_aln.fasta","w") as all_DNA_seqs:
        with open(output_folder / "all_AA_aln.fasta","w") as all_AA_seqs:
            for index, row in df.iterrows():
                seqid = row['Id']
                seqname = row['Name']
                all_DNA_seqs.write(">{}\n".format(seqname))
                all_DNA_seqs.write("{}\n".format(DNA_sequences[str(seqid)].seq))
                all_AA_seqs.write(">{}\n".format(seqname))
                all_AA_seqs.write("{}\n".format(AA_sequences[str(seqid)].seq))

def analyse_DNA_sequence(sequence,df,weight=1):
    '''
    Read the DNA codons of the given sequence and add them to the dataframe with the given ponderation.
    '''
    for i in range(len(sequence)//3):
        codon = sequence[i*3:(i+1)*3]
        if(codon in Codons_list):
            df.at[i,codon] += weight
        else:
            df.at[i,"Other"] += weight
    return df

def analyse_AA_sequence(sequence,df,weight=1):
    '''
    Read the AA residues of the given sequence and add them to the dataframe with the given ponderation.
    '''
    for i in range(len(sequence)):
        aa = sequence[i]
        if(aa in AA_list):
            df.at[i,aa] += weight
        else:
            df.at[i,"Other"] += weight
    return df

def analyse_alignment(output_folder):
    '''
    Read DNA and AA alignments and summarize them in four dataframes (AA+DNA sequences / Unique+All sequences).
    '''
    df = pd.read_csv(output_folder / "sequences" / "sequence_unique_id.csv", sep=";")
    df = df.groupby("Id").count().reset_index()
    df.columns = ["Id","Count"]
    DNA_sequences = fasta2dict(output_folder / "sequences" / "unique_DNA_aln.fasta")
    AA_sequences = fasta2dict(output_folder / "sequences" / "unique_AA_aln.fasta")
    first_dna_seq = list(DNA_sequences.values())[0]
    n_DNA = len(first_dna_seq)
    first_aa_seq = list(AA_sequences.values())[0]
    n_AA = len(first_aa_seq)
    unique_DNA_df = pd.DataFrame(np.zeros((n_DNA//3,len(Codons_list))))
    unique_DNA_df.columns = Codons_list
    unique_AA_df = pd.DataFrame(np.zeros((n_AA,len(AA_list))))
    unique_AA_df.columns = AA_list
    all_DNA_df = pd.DataFrame(np.zeros((n_DNA//3,len(Codons_list))))
    all_DNA_df.columns = Codons_list
    all_AA_df = pd.DataFrame(np.zeros((n_AA,len(AA_list))))
    all_AA_df.columns = AA_list
    for index, row in df.iterrows():
        seqid = str(row['Id'])
        seqcount = row['Count']
        DNA_seq = str(DNA_sequences[seqid].seq).upper()
        AA_seq = str(AA_sequences[seqid].seq).upper()
        unique_DNA_df = analyse_DNA_sequence(DNA_seq,unique_DNA_df)
        all_DNA_df = analyse_DNA_sequence(DNA_seq,all_DNA_df,weight=seqcount)
        unique_AA_df = analyse_AA_sequence(AA_seq,unique_AA_df)
        all_AA_df = analyse_AA_sequence(AA_seq,all_AA_df,weight=seqcount)
    unique_DNA_df.to_csv(output_folder / "unique_DNA.csv",sep=";",index=None)
    all_DNA_df.to_csv(output_folder / "all_DNA.csv" ,sep=";",index=None)
    unique_AA_df.to_csv(output_folder / "unique_AA.csv",sep=";",index=None)
    all_AA_df.to_csv(output_folder / "all_AA.csv",sep=";",index=None)

def homologs_analysis(output_folder):
    '''
    Entire workflow for parallelization.
    '''
    extract_unique_sequences(output_folder/ "homologs_analysis", output_folder/ "reference_analysis" / "reference_sequence.fasta")
    translate_unique_sequences(output_folder / "homologs_analysis" / "sequences")
    align_AA_sequences(output_folder / "homologs_analysis" / "sequences")
    align_DNA_sequences(output_folder / "homologs_analysis" / "sequences")
    remove_reference(output_folder / "homologs_analysis" / "sequences" / "unique_AA_aln.fasta")
    remove_reference(output_folder / "homologs_analysis" / "sequences" / "unique_DNA_aln.fasta")
    retrieve_all_sequences(output_folder / "homologs_analysis" / "sequences")
    analyse_alignment(output_folder / "homologs_analysis")

def distant_homologs_analysis(output_folder):
    '''
    Read AA alignment and summarize it in a dataframe, compute CI entropy as well.
    '''
    AA_sequences = read_fasta(output_folder / "distant_homologs.fasta")
    n_AA = len(AA_sequences[0])
    AA_df = pd.DataFrame(np.zeros((n_AA,len(AA_list))))
    AA_df.columns = AA_list
    for AA_seq in AA_sequences:
        AA_df = analyse_AA_sequence(str(AA_seq.seq).upper(),AA_df)
    AA_df.to_csv(output_folder / "distant_homologs.csv",sep=";",index=None)
    AA_df["Sum"] = AA_df[AA_list[:-2]].apply(sum,axis=1)
    for character in AA_list[:-2]:
        AA_df["{}_proba".format(character)] = AA_df[character]/AA_df["Sum"]
        AA_df["{}_entropy".format(character)] = AA_df["{}_proba".format(character)].apply(entropy)
    entropy_columns = ["{}_entropy".format(character) for character in AA_list[:-2]]
    AA_df["CI_Entropy"] = AA_df[entropy_columns].apply(sum,axis=1)
    AA_df[["CI_Entropy"]].to_csv(output_folder / "CI_Entropy.csv",sep=";",index=None)

def gather_possible_mutants(directory):
    '''
    For each gene, gather scores of possible mutants + information on the type of mutation (accessible, profile).
    '''
    gene = directory.stem
    profile = pd.read_csv(directory / "distant_homologs" / "distant_homologs.csv", sep=";")
    accessible = pd.read_csv(directory / "reference_analysis" / "accessible_mutants.csv")
    observed = pd.read_csv(directory / "homologs_analysis" / "all_AA.csv",sep=";")
    DCA = pd.read_csv(directory / "reference_analysis" / "DCA_scores.csv")
    codons = list(accessible.columns)
    codons.remove('Codon')
    codons.remove('*')
    observed = observed.loc[:,codons]
    obs_total = observed.sum(axis=1)
    accessible = accessible.loc[:,codons]>0
    profile = profile.loc[:,codons]
    profile_total = profile.sum(axis=1)
    for aa in codons:
        observed[aa] = observed[aa]/obs_total
        profile[aa] = profile[aa]/profile_total
    DCA_aa = DCA.loc[:,codons]
    Locus, AA, Ref, Observed, Profile, Accessible, Score = list(), list(), list(), list(), list(), list(), list()
    for i in list(DCA_aa.stack().index):
        reference_AA = DCA.loc[(i[0],"Context")]
        Ref.append(i[1]==reference_AA)
        Locus.append(i[0])
        AA.append(i[1])
        Observed.append(observed.loc[(i[0],i[1])])
        Profile.append(profile.loc[(i[0],i[1])])
        Accessible.append(accessible.loc[(i[0],i[1])])
        Score.append(DCA.loc[(i[0],i[1])])
    results = pd.DataFrame({"Locus":Locus,"AA":AA,"Reference":Ref,"Observed":Observed,"Profile":Profile,"Accessible":Accessible,"DCA_score":Score})
    results["Gene"] = gene
    results[["Gene","Locus","Reference","AA","Observed","Profile","Accessible","DCA_score"]].to_csv(directory / "stats" / "mutants_sites.csv", index=None)

def site_freq(directory):
    '''
    Write csv with stats on polymorphism frequencies in homologs and distant homologs (gaps, others, maf, reference allele freq).
    '''
    gene = (directory.stem)
    DCA_proba = pd.read_csv(directory / "reference_analysis" / "DCA_proba.csv",sep=",")
    distant_homologs = pd.read_csv(directory / "distant_homologs" / "distant_homologs.csv",sep=";")
    all_AA = pd.read_csv(directory / "homologs_analysis" / "all_AA.csv",sep=";")
    distant_homologs["Total"] = distant_homologs.sum(axis=1)
    distant_homologs["Prop_gaps_distant"] = (distant_homologs["-"]+distant_homologs["Other"])/distant_homologs["Total"]
    distant_homologs["Gene"] = gene
    distant_homologs["Locus"] = -1
    for i in range(len(distant_homologs)):
        distant_homologs.at[i,"Locus"] = i
    col_list = list(all_AA)
    col_list.remove('-')
    col_list.remove('Other')
    all_AA["Total"] = all_AA.sum(axis=1)
    all_AA["Norm_total"] = all_AA[col_list].sum(axis=1)
    all_AA["Max"] = all_AA[col_list].max(axis=1)
    all_AA["Prop_gaps_local"] = (all_AA["-"]+all_AA["Other"])/all_AA["Total"]
    all_AA["Freq"] = 1-all_AA["Max"]/all_AA["Norm_total"]
    all_AA["Gene"] = gene
    all_AA["Locus"] = -1
    all_AA["Context_freq"] = -1
    for i in range(len(all_AA)):
        aa = DCA_proba["Context"][i]
        all_AA.at[i,"Locus"] = i
        all_AA.at[i,"Context_freq"] = all_AA[aa][i]
    all_AA["Context_freq"] = all_AA["Context_freq"]/all_AA["Total"]
    all_AA[["Gene","Locus","Prop_gaps_local","Freq","Context_freq"]].merge(distant_homologs[["Gene","Locus","Prop_gaps_distant"]]).to_csv(directory / "stats" / "site_frequencies.csv", index=None)

def stats_entropy(directory):
    '''
    Gather CD and CI entropy in same csv.
    '''
    gene = (directory.stem)
    CI_entropy = pd.read_csv(directory / "distant_homologs" / "CI_Entropy.csv",sep=",")
    CD_entropy = pd.read_csv(directory / "reference_analysis" / "CD_Entropy.csv",sep=",")
    CI_entropy["Gene"] = gene
    CI_entropy["Locus"] = -1
    for i in range(len(CI_entropy)):
        CI_entropy.at[i,"Locus"] = i
    CD_entropy["Gene"] = gene
    CD_entropy["Locus"] = -1
    for i in range(len(CD_entropy)):
        CD_entropy.at[i,"Locus"] = i
    CD_entropy[["Gene","Locus","CD_Entropy"]].merge(CI_entropy[["Gene","Locus","CI_Entropy"]]).to_csv(directory / "stats" / "entropies.csv", index=None)

def stats_proba_polymorphisms(directory, threshold):
    '''
    DCA proba of observed polymorphisms at frequencies over the given threshold.
    '''
    gene = (directory.stem)
    DCA_proba = pd.read_csv(directory / "reference_analysis" / "DCA_proba.csv",sep=",")
    all_AA = pd.read_csv(directory / "homologs_analysis" / "all_AA.csv",sep=";")
    accessible = pd.read_csv(directory / "reference_analysis" / "accessible_mutants.csv",sep=",")
    AA_col = list(all_AA)
    AA_col.remove("-")
    AA_col.remove("Other")
    DCA_col = list(DCA_proba)
    DCA_col.remove("Context")
    accessible = accessible[AA_col].applymap(lambda x: int(x>0))
    accessible = accessible.loc[0:len(DCA_proba),:]
    all_AA = all_AA[AA_col]
    all_AA["Threshold"] = all_AA.sum(axis=1)*threshold
    for col in AA_col:
        all_AA[col] = (all_AA[col]>all_AA["Threshold"]).astype('int')
    all_AA["Proba_obs"] = np.multiply(DCA_proba[DCA_col],all_AA[AA_col]).sum(axis=1)
    all_AA["Proba_obs_acc"] = np.multiply(np.multiply(DCA_proba[DCA_col],all_AA[AA_col]),accessible[AA_col]).sum(axis=1)
    all_AA["Tot_acc"] = (accessible[AA_col]).sum(axis=1)
    all_AA["Proba_tot_acc"] = np.multiply(DCA_proba[DCA_col],accessible[AA_col]).sum(axis=1)
    all_AA["Gene"] = gene
    all_AA["Locus"] = -1
    for i in range(len(DCA_proba)):
        all_AA.at[i,"Locus"] = i
    filename = "proba_polymorphisms_{}.csv".format(threshold)
    all_AA[["Gene","Locus","Proba_obs","Proba_obs_acc","Proba_tot_acc","Tot_acc"]].to_csv(directory / "stats" / filename, index=None)

def stats_synonymous_polymorphisms(directory):
    '''
    Number of synonymous mutants per site (observed, accessible,...).
    '''
    gene = (directory.stem)
    DCA_proba = pd.read_csv(directory / "reference_analysis" / "DCA_proba.csv",sep=",")
    accessible = pd.read_csv(directory / "reference_analysis" / "accessible_codons.csv",sep=",")
    synonymous = pd.read_csv(directory / "reference_analysis" / "synonymous_mutants.csv",sep=",")
    all_DNA = pd.read_csv(directory / "homologs_analysis" / "all_DNA.csv",sep=";")
    codons = list(synonymous)
    codons.remove("Codon")
    all_DNA = all_DNA[codons].applymap(lambda x: int(x>0))
    synonymous["Locus"] = -1
    all_DNA["Locus"] = -1
    accessible["Locus"] = -1
    for i in range(len(DCA_proba["Context"])):
        synonymous.at[i,"Locus"] = i
        all_DNA.at[i,"Locus"] = i
        accessible.at[i,"Locus"] = i
    all_DNA = all_DNA[all_DNA["Locus"]>-1]
    synonymous = synonymous[synonymous["Locus"]>-1]
    accessible = accessible[accessible["Locus"]>-1] 
    synonymous["Observed_syn"] = (np.multiply(all_DNA[codons],synonymous[codons]).sum(axis=1)).astype('int')
    synonymous["Observed_acc_syn"] = (np.multiply(accessible[codons],np.multiply(all_DNA[codons],synonymous[codons])).sum(axis=1)).astype('int')
    synonymous["Total_syn"] = (synonymous[codons].sum(axis=1)).astype('int')
    synonymous["Total_acc_syn"] = (np.multiply(accessible[codons],synonymous[codons]).sum(axis=1)).astype('int')
    synonymous["Observed_nsyn"] = (np.multiply(all_DNA[codons],1-synonymous[codons]).sum(axis=1)).astype('int')
    synonymous["Observed_acc_nsyn"] = (np.multiply(accessible[codons],np.multiply(all_DNA[codons],1-synonymous[codons])).sum(axis=1)).astype('int')
    synonymous["Total_nsyn"] = ((1-synonymous[codons]).sum(axis=1)).astype('int')
    synonymous["Total_acc_nsyn"] = (np.multiply(accessible[codons],1-synonymous[codons]).sum(axis=1)).astype('int')
    synonymous["Gene"] = gene
    synonymous[["Gene","Locus","Observed_syn","Observed_acc_syn","Total_syn","Total_acc_syn","Observed_nsyn","Observed_acc_nsyn","Total_nsyn","Total_acc_nsyn"]].to_csv(directory / "stats" / "nsynonymous_synonymous_polymorphisms.csv", index=None)

def stats_count_polymorphisms(directory, threshold):
    '''
    DCA proba of observed polymorphisms at frequencies over the given threshold.
    '''
    gene = (directory.stem)
    DCA_proba = pd.read_csv(directory / "reference_analysis" / "DCA_proba.csv",sep=",")
    all_AA = pd.read_csv(directory / "homologs_analysis" / "all_AA.csv",sep=";")
    accessible = pd.read_csv(directory / "reference_analysis" / "accessible_mutants.csv",sep=",")
    AA_col = list(all_AA)
    AA_col.remove("-")
    AA_col.remove("Other")
    DCA_col = list(DCA_proba)
    DCA_col.remove("Context")
    Reference_AA = DCA_proba["Context"]
    accessible = accessible[AA_col].applymap(lambda x: int(x>0))
    accessible = accessible.loc[0:len(DCA_proba),:]
    all_AA = all_AA[AA_col]
    for index, row in DCA_proba.iterrows():
        col = "{}_proba".format(Reference_AA[index])
        DCA_proba.loc[index,col] = 0
    DCA_proba = DCA_proba[DCA_col].applymap(lambda x: int(x>threshold))
    all_AA = all_AA[AA_col].applymap(lambda x: int(x>0))
    all_AA["Count_obs"] = np.multiply(DCA_proba[DCA_col],all_AA[AA_col]).sum(axis=1)
    all_AA["Count_obs_acc"] = np.multiply(np.multiply(DCA_proba[DCA_col],all_AA[AA_col]),accessible[AA_col]).sum(axis=1)
    all_AA["Tot_dca"] = (DCA_proba[DCA_col]).sum(axis=1)
    all_AA["Tot_acc_dca"] = np.multiply(DCA_proba[DCA_col],accessible[AA_col]).sum(axis=1)
    all_AA["Gene"] = gene
    all_AA["Locus"] = -1
    for i in range(len(DCA_proba)):
        all_AA.at[i,"Locus"] = i
    filename = "count_polymorphisms_{}.csv".format(threshold)
    all_AA[["Gene","Locus","Count_obs","Count_obs_acc","Tot_dca","Tot_acc_dca"]].to_csv(directory / "stats" / filename, index=None)

def compute_site_stats(directory, threshold):
    '''
    Compute all stats for the given gene directory, using the specified threshold for polymorphism.
    '''
    site_freq(directory)
    stats_entropy(directory)
    stats_proba_polymorphisms(directory, threshold)
    stats_synonymous_polymorphisms(directory)
    stats_count_polymorphisms(directory, 0.05)

def merge_dfs(directory, threshold):
    '''
    Merge all stats from each individual gene into one dataframe.
    '''
    filename = "proba_polymorphisms_{}.csv".format(threshold)
    proba_polymorphisms = pd.read_csv(directory / "stats" / filename)
    count_polymorphisms = pd.read_csv(directory / "stats" / "count_polymorphisms_0.05.csv")
    synonymous = pd.read_csv(directory / "stats" / "nsynonymous_synonymous_polymorphisms.csv")
    entropies = pd.read_csv(directory / "stats" / "entropies.csv")
    site_freq = pd.read_csv(directory / "stats" / "site_frequencies.csv")
    return site_freq.merge(entropies, on=["Gene","Locus"], how="left").merge(proba_polymorphisms, on=["Gene","Locus"], how="left").merge(synonymous, on=["Gene","Locus"], how="left").merge(count_polymorphisms, on=["Gene","Locus"], how="left")

def gather_info_sites(directory):
    '''
    For each gene, create a dataframe for JC69 simulations.
    '''
    gene = (directory.stem)
    accessible = pd.read_csv(directory / "reference_analysis" / "accessible_mutants.csv")
    DCA_proba = pd.read_csv(directory / "reference_analysis" / "DCA_proba.csv")
    codons = list(accessible.columns)
    codons.remove('Codon')
    codons.remove('*')
    df = accessible.loc[:,codons]
    Locus = list()
    AA = list()
    Counts = list()
    Acceptance_th = list()
    Proba = list()
    Ref = list()
    for i in list(df[df>0].stack().index):
        reference_AA = DCA_proba.loc[(i[0],"Context")]
        is_ref = int(i[1]==reference_AA)
        Locus.append(i[0]+1)
        AA.append(i[1])
        Counts.append(accessible.loc[i]-is_ref)
        reference_proba = DCA_proba.loc[(i[0],"{}_proba".format(reference_AA))]
        AA_proba = DCA_proba.loc[(i[0],"{}_proba".format(i[1]))]
        Acceptance_th.append(AA_proba/(AA_proba+reference_proba))
        Proba.append(AA_proba)
        Ref.append(reference_AA)
    results = pd.DataFrame({"Locus":Locus,"AA":AA,"Counts":Counts,"Acceptance_th":Acceptance_th,"Proba":Proba,"Ref":Ref})
    results["Gene"] = gene
    results[["Gene","Locus","AA","Counts","Acceptance_th","Proba","Ref"]].to_csv(directory / "stats" / "simulated_sites.csv", index=None)



def single_muts(seq,refseq):
    '''
    Compare seq to refseq and look for amino acid differences.
    The input sequences should have exactly 2 differences, if these are not gaps,
    the function returns the corresponding single mutants (2 sequences with 1 amino
    acid change each) and the loci of the fixed differences.
    '''
    loci = [i for i in range(len(refseq)) if seq[i]!=refseq[i]]
    single_mut_1, single_mut_2 = list(refseq), list(refseq)
    aa1, aa2 = seq[loci[0]],  seq[loci[1]]
    if(aa1 != "-" and aa2 != "-"):
        single_mut_1[loci[0]] = aa1
        single_mut_2[loci[1]] = aa2
        return ["".join(single_mut_1), "".join(single_mut_2), loci[0], loci[1]]
    else:
        return []

def double_muts(dca_models_path, folder):
    '''
    Compute score of double mutants and compare it to the sum of the scores of single mutants
    for sequences that have exactly 2 amino acid changes from the reference.
    '''
    refseq = translate(read_fasta(folder / "reference_analysis" / "reference_sequence.fasta")[0].seq)
    msa = [str(sequence.seq) for sequence in read_fasta(folder / "homologs_analysis" / "sequences" / "unique_AA_aln.fasta")]
    selection = [seq for seq in msa if sum([int(seq[i]!=refseq[i]) for i in range(len(refseq))])==2]
    sequences = {seq:single_muts(seq,refseq) for seq in selection}
    h,J = read_dca_par(dca_models_path, folder.stem)
    E_refseq = compute_energy(refseq, h, J)
    Double_mut, Sum_single_mut, Loci_1, Loci_2 = list(), list(), list(), list()
    if(all_standard_aa(refseq)):
        for key in sequences.keys():
            values = sequences[key]
            if(len(values)==4):
                seq_1, seq_2, loci_1, loci_2 = values[0], values[1], values[2], values[3]
                if(all_standard_aa(key)):
                    dble_mut = compute_energy(key, h, J)-E_refseq
                    sgle_mut = compute_energy(seq_1, h, J)+compute_energy(seq_2, h, J)-2*E_refseq
                    Double_mut.append(dble_mut)
                    Sum_single_mut.append(sgle_mut)
                    Loci_1.append(loci_1)
                    Loci_2.append(loci_2)
    pd.DataFrame({"Double_mut":Double_mut,"Sum_single_mut":Sum_single_mut, "Loci_1":Loci_1, "Loci_2":Loci_2}).to_csv(folder / "stats" / "double_mut_epistasis.csv", index=None)

def merge_df(filename, input_folder, output_folder):
    '''
    Merge all individual gene dataframes written in the different
    `filename` files (one per gene) into one common file.
    '''
    dfs = list()
    for folder in input_folder.iterdir():
        if((folder / "stats" / filename).is_file()):
            df = pd.read_csv(folder / "stats" / filename)
            df["Gene"] = folder.stem
            dfs.append(df)
    pd.concat(dfs).to_csv(output_folder / filename, index=None)

def IPR(dca_models_path, folder):
    '''
    Compute IPR values for all positions in the reference sequence corresponding to the input folder.
    '''
    refseq = translate(read_fasta(folder / "reference_analysis" / "reference_sequence.fasta")[0].seq)
    h,J = read_dca_par(dca_models_path, folder.stem)
    Locus, IPR = list(), list()
    for i in range(len(refseq)):
        J_vector = list()
        for j in range(i):
            J_vector.append(J[d_aa_num[refseq[j]], d_aa_num[refseq[i]], j, i])
        for j in range(i+1, len(refseq)):
            J_vector.append(J[d_aa_num[refseq[i]], d_aa_num[refseq[j]], i, j])
        J_vector = np.array(J_vector)
        Locus.append(i)
        IPR.append(np.sum((J_vector**2/np.sum(J_vector**2))**2))
    pd.DataFrame({"Locus":Locus, "IPR":IPR}).to_csv(folder / "stats" / "IPR.csv", index=None)






def aa_to_df(sequence):
    '''
    Change a AA sequence to a dataframe with 1 for AA in the sequence and 0 elsewhere.
    '''
    seq_up = sequence.upper()
    AAs = valid_aa + ['*','X'] # MODIFIED
    n = len(seq_up)
    df = pd.DataFrame(np.zeros((n,len(AAs))),columns=AAs)
    for i in range(n):
        df.loc[i,seq_up[i]] += 1
    return df

def dcamatrix2df(matrix, sequence, output_folder, filename):
    '''
    Compute quantities from DCA scores and write results in corresponding csv files.
    '''
    df = pd.DataFrame(matrix)
    df.columns = valid_aa
    for character in valid_aa[:-1]:
        df["{}_exp".format(character)] = np.exp(-df[character])
    exp_columns = ["{}_exp".format(character) for character in valid_aa[:-1]]
    df["Sum"] = df[exp_columns].apply(sum,axis=1)
    for character in valid_aa[:-1]:
        df["{}_proba".format(character)] = df["{}_exp".format(character)]/df["Sum"]
        df["{}_entropy".format(character)] = df["{}_proba".format(character)].apply(entropy)
    proba_columns = ["{}_proba".format(character) for character in valid_aa[:-1]]
    entropy_columns = ["{}_entropy".format(character) for character in valid_aa[:-1]]
    df["CD_Entropy"] = df[entropy_columns].apply(sum,axis=1) 
    df[["CD_Entropy"]].to_csv(output_folder / "CD_Entropy" / filename, index=None)
    df[valid_aa].to_csv(output_folder / "DCA_scores" / filename, index=None)
    df[proba_columns].to_csv(output_folder / "DCA_proba" / filename, index=None, header=valid_aa[:-1])

def distant_homologs2dfs(AA_sequences, output_folder):
    '''
    Read AA alignment and summarize it in a dataframe, compute CI entropy as well.
    '''
    AA_list = valid_aa + ['Other']
    n_AA = len(AA_sequences[0])
    AA_df = pd.DataFrame(np.zeros((n_AA,len(AA_list))))
    AA_df.columns = AA_list
    for AA_seq in AA_sequences:
        AA_df = analyse_AA_sequence(str(AA_seq.seq).upper(),AA_df)
    AA_df.to_csv(output_folder / "Polymorphism" / "distant_homologs.csv",index=None)
    AA_df["Sum"] = AA_df[AA_list[:-2]].apply(sum,axis=1)
    for character in AA_list[:-2]:
        AA_df["{}_proba".format(character)] = AA_df[character]/AA_df["Sum"]
        AA_df["{}_entropy".format(character)] = AA_df["{}_proba".format(character)].apply(entropy)
    entropy_columns = ["{}_entropy".format(character) for character in AA_list[:-2]]
    AA_df["CI_Entropy"] = AA_df[entropy_columns].apply(sum,axis=1)
    AA_df[["CI_Entropy"]].to_csv(output_folder / "CI_Entropy.csv",index=None)

def make_output_folders(gene, output_folder):
    '''
    Create a gene folder with subfolders in the given output folder to store
    temporary files.
    '''
    (output_folder / gene / "CD_Entropy").mkdir(parents=True, exist_ok=True)
    (output_folder / gene / "DCA_scores").mkdir(parents=True, exist_ok=True)
    (output_folder / gene / "DCA_proba").mkdir(parents=True, exist_ok=True)
    (output_folder / gene / "Polymorphism").mkdir(parents=True, exist_ok=True)

def process_gene(gene, output_folder, msa_folder, aln_folder, dca_path):
    '''
    Read sequences of closely related species, store fixed differences and compute DCA scores.
    '''
    make_output_folders(gene, output_folder)
    aln_file = aln_folder / f"{gene}.fa"
    msa_file = msa_folder / f"{gene}.fa"
    h,J = read_dca_par(dca_path, gene)
    distant_sequences = read_fasta(aln_file)
    seqid = list()
    energy = list()
    for sequence in distant_sequences:
        filename = "{}.csv".format(sequence.id)
        AA_seq = str(sequence.seq)
        aa_to_df(AA_seq).to_csv(output_folder / gene / "Polymorphism" / filename, index=None)
        energy_matrix = perform_DCA(AA_seq, J, h)
        dcamatrix2df(energy_matrix, AA_seq, output_folder / gene, filename)
        energy.append(compute_energy(AA_seq, h, J))
        seqid.append(sequence.id)
    pd.DataFrame({"Sequence_id":seqid,"Energy":energy}).to_csv(output_folder / gene / "sequence_scores.csv",index=None)
    distant_homologs = read_fasta(msa_file)
    distant_homologs2dfs(distant_homologs, output_folder / gene)
            
def full_seq(folder):
    '''
    Compute total epistatic cost for pairs of homologuous sequences of closely 
    related species (if they do not have more than 1-gap difference).
    '''
    AAs = [aa for aa in valid_aa if aa!="-"]
    full_seq_dca = pd.read_csv(folder / "sequence_scores.csv")
    full_seq_dca = full_seq_dca.set_index("Sequence_id")
    Species_i, Species_j, FD, Single_muts, Full_seq = list(), list(), list(), list(), list()
    species_files = [file.name for file in (folder / "DCA_scores").iterdir() if file.suffix==".csv"]
    for i in range(len(species_files)):
        for j in range(i+1,len(species_files)):
            polymorphism_i = pd.read_csv(folder / "Polymorphism" / species_files[i])
            polymorphism_j = pd.read_csv(folder / "Polymorphism" / species_files[j])
            gaps = np.multiply(polymorphism_i["-"],1-polymorphism_j["-"]).sum() + np.multiply(1-polymorphism_i["-"], polymorphism_j["-"]).sum()
            if(gaps <= 1):
                dca_i = pd.read_csv(folder / "DCA_scores" / species_files[i])
                dca_j = pd.read_csv(folder / "DCA_scores" / species_files[j])
                fixed_differences = np.multiply(polymorphism_i[AAs],1-polymorphism_j[AAs]).sum().sum()/len(polymorphism_i)
                sum_single_muts_i = np.multiply(polymorphism_i[AAs+["-"]],dca_j[AAs+["-"]]).sum().sum()
                sum_single_muts_j = np.multiply(polymorphism_j[AAs+["-"]],dca_i[AAs+["-"]]).sum().sum()
                species_i, species_j = species_files[i].split(".")[0], species_files[j].split(".")[0]
                Species_i.append(species_i)
                Species_j.append(species_j)
                energy = full_seq_dca.loc[species_j,"Energy"]-full_seq_dca.loc[species_i,"Energy"]
                FD.append(fixed_differences)
                Single_muts.append(sum_single_muts_i)
                Full_seq.append(-energy)
                Species_i.append(species_files[j].split(".")[0])
                Species_j.append(species_files[i].split(".")[0])
                FD.append(fixed_differences)
                Single_muts.append(sum_single_muts_j)
                Full_seq.append(energy)
    df = pd.DataFrame({"Species_i":Species_i,"Species_j":Species_j,"FD":FD,"Single_muts":Single_muts,"Full_seq":Full_seq})
    df.to_csv(folder/"full_seq_single_muts.csv",index=None)

def merge_full_seq_dfs(data_folder, output_file):
    '''
    Merge all individual gene dataframes written in the different
    full_seq_single_muts.csv files (one per gene) into one common file.
    '''
    dfs = list()
    for folder in data_folder.iterdir():
        if((folder/"full_seq_single_muts.csv").is_file()):
            df = pd.read_csv(folder/"full_seq_single_muts.csv")
            df["Gene"] = folder.stem
            dfs.append(df)
    pd.concat(dfs).to_csv(output_file, index=None)

def compare_sequences(sequence_1, sequence_2):
    '''
    Retrieve all pairs of fixed differences between 2 sequences.
    '''
    loci = [i for i in range(len(sequence_1)) if sequence_1[i]!=sequence_2[i]]
    I, K, AAi1, AAi2, AAk1, AAk2= list(), list(), list(), list(), list(), list()
    for i in range(len(loci)):
        aai1 = sequence_1[loci[i]]
        aai2 = sequence_2[loci[i]]
        for k in range(i+1,len(loci)):
            aak1 = sequence_1[loci[k]]
            aak2 = sequence_2[loci[k]]
            I.append(loci[i])
            K.append(loci[k])
            AAi1.append(aai1)
            AAi2.append(aai2)
            AAk1.append(aak1)
            AAk2.append(aak2)
    return pd.DataFrame({"i":I, "k":K, "AAi1":AAi1, "AAi2":AAi2, "AAk1":AAk1, "AAk2":AAk2})

def retrieve_couplings(df, J, seq): #MODIFIED
    '''
    Find DCA epistatic couplings between pairs of amino acids in df.
    '''
    df["J_i1k1"] = None
    df["J_i2k2"] = None
    df["J_i1k2"] = None
    df["J_i2k1"] = None
    for index, row in df.iterrows():
        df.at[index,"J_i1k1"] = -J[d_aa_num[row["AAi1"]], d_aa_num[row["AAk1"]], row["i"], row["k"]]
        df.at[index,"J_i2k2"] = -J[d_aa_num[row["AAi2"]], d_aa_num[row["AAk2"]], row["i"], row["k"]]
        df.at[index,"J_i1k2"] = -J[d_aa_num[row["AAi1"]], d_aa_num[row["AAk2"]], row["i"], row["k"]]
        df.at[index,"J_i2k1"] = -J[d_aa_num[row["AAi2"]], d_aa_num[row["AAk1"]], row["i"], row["k"]]
    return df

def couplings(gene, dca_path, aln_folder):
    '''
    Retrieve epistatic couplings between amino acid changes
    fixed between E. coli and Y. pestis.
    '''
    h,J = read_dca_par(dca_path, gene.stem)
    fasta_filename = "{}.fa".format(gene.stem)
    distant_sequences = read_fasta(aln_folder / fasta_filename)
    for sequence in distant_sequences:
        if(sequence.id=="ESC_GA4805AA"):
            refseq = str(sequence.seq)
        elif(sequence.id=="Yersiniia_pestis"):
            pestis = str(sequence.seq)
    try:
        retrieve_couplings(compare_sequences(refseq, pestis), J, refseq).to_csv(gene / "couplings.csv", index=None)
    except:
        print("Error: {}".format(gene.stem))

def merge_couplings_dfs(data_folder, output_file):
    '''
    Merge all individual gene dataframes written in the different
    couplings.csv files (one per gene) into one common file.
    '''
    dfs = list()
    for folder in data_folder.iterdir():
        if((folder / "couplings.csv").is_file()):
            df = pd.read_csv(folder / "couplings.csv")
            df["Gene"] = folder.stem
            dfs.append(df)
    pd.concat(dfs).to_csv(output_file, index=None)
