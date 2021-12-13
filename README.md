# DCA to decipher polymorphism in E. coli strains
We use computational models based on Direct Coupling Analysis - [DCA](https://en.wikipedia.org/wiki/Direct_coupling_analysis) - trained on [PFAM](http://pfam.xfam.org/) domains of distant distant homologues to accurately predict the polymorphisms segregating in a panel of 61,157 Escherichia coli genomes. 

We show that the genetic context (i.e. the rest of the protein sequence) strongly constrains the tolerable amino acids in 30% to 50% of amino-acid sites. Our study also suggests the gradual build-up of genetic context over long evolutionary timescales by the accumulation of small epistatic contributions.


Paper: [Deciphering polymorphism in 61,157 Escherichia coli genomes via epistatic sequence landscapes](link to the paper) (Vigué L.\*,  Croce G.\*, and al. xxxxxxx , 2021)

![figure](ecoli_sequence_landscape.png)

We provide here the code to reproduce the key results and figures of the paper.

## Installation:
To run the code, you first need to install :
- python3: (the code was tested on python v3.8)
- [julia:](https://julialang.org/)  to run the DCA pseudo-likelihood inference algorithm  (tested on julia v1.6)
- [mafft:](https://mafft.cbrc.jp/alignment/software/) to align sequences (tested on v7.471 (2020/Jul/3))

Then clone the repository to a directory of your choice, where you have writing permissions, and install the python libraries by running:
```
pip install requirements.txt
```
It is strongly recommended to use a virtual environment.
 
You also need to install plmDCA (pseudo-likelihood inference algorithm) for julia (see how to do it from https://github.com/pagnani/PlmDCA)

## Config your paths:
Open the file  ```src/config.py``` with your favorite editor, and replace ```path_julia``` with the path to the julia executable on your computer.

## Usage:
Our aim is to train a DCA model on distant homologues (PFAM data - long term evolution - highly variable sequences varibility) and use it to predict polymorphism in E. coli strains (short term evolution - most positions are highly conserved).

## Reproduce key results

Download data from Zenodo at https://zenodo.org/record/5774192#.YbUZILvjLJE (DOI 10.5281/zenodo.5774191) and put the tar archive in this repository. Then use following commands to perform data analysis.

```
./extract_datasets.sh
python3 train_dca_models.py
python3 analyse_coli_strains.py
python3 analyse_closely_diverged_species.py
jupyter lab Produce_Figures.ipynb
```
