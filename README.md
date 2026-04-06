# Evolutionary of histone mRNA 3' end processing machinery in eukaryotes

A comprehensive pipeline for identifying and analyzing Stem-Loop Binding Protein (SLBP) orthologs and replication-dependent histone processing elements across eukaryotes.

## Repository Structure
Histone_SLBP_Project/
|- bin/	#Executable Bash and Python Scripts
|- data/	#Raw input FASTA files
|- results/	#Analysis outputs, e.g. metafile, figures and etc. 
|- logs/	#errors
|- requirements.txt	#Packages required for this project
|- Readme.md

## Pipeline Overview
1. **Protein Identification**: Remote BLASTP and jackhmmer searches for SLBP, LSM10/11, and FLASH.
2. **Domain Verification**: RPS-BLAST and HMMER profiling for RBD (PF15247) and HFD domains.
3. **RNA Analysis**: Detection of histone 3' UTR Stem-Loops and PAS motifs.
4. **Evolutionary Synthesis**: Gene-species tree reconciliation and processing state classification.

## Requirements
- NCBI BLAST+
- HMMER 3.x
- RNAfold (ViennaRNA)
- Python 3.8+
- cmsearch
- MAFFT
- FastTree
- ete3


