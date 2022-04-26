using DelimitedFiles
using DCAUtils
using Distributed

filelist = readdir("./datasets/DCA_training_MSAs/local_strains/")
mkdir("./datasets/weighted_frequencies")


#### settings
max_gap_fraction = 0.9
remove_dups = true
θ = 0.2

#@distributed for filename in filelist
for filename in filelist
   if endswith(filename, ".fa")
       aln = DCAUtils.ReadFastaAlignment.read_fasta_alignment("./datasets/DCA_training_MSAs/local_strains/"*filename, max_gap_fraction)
       if remove_dups
           aln, _ = remove_duplicate_sequences(aln)
           N, M = size(aln)
           q = Int(maximum(aln))
           Pi_true, Pij_true, Meff, _ = compute_weighted_frequencies(aln, q, θ)
       end
       open("./datasets/weighted_frequencies/"*filename, "w") do io
           writedlm(io, ["A" "C" "D" "E" "F" "G" "H" "I" "K" "L" "M" "N" "P" "Q" "R" "S" "T" "V" "W" "Y"], ",")
           writedlm(io, transpose(reshape(Pi_true, (20,N))), ",") 
       end
   end
end

