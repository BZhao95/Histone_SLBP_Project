# SLBP and Histone Processing

A comprehensive pipeline for identifying Stem-Loop Binding Protein (SLBP) orthologs and replication-dependent histone processing elements across eukaryotes.

## Pipeline Overview
1. **Protein Identification**: Remote BLASTP and jackhmmer searches for SLBP, LSM10/11, and FLASH.
2. **Domain Verification**: RPS-BLAST and HMMER profiling for RBD (PF15247) and HFD domains.
3. **RNA Analysis**: Detection of histone 3' UTR Stem-Loops and PAS motifs.
4. **Evolutionary Synthesis**: Gene-species tree reconciliation and processing state classification.

## Requirements
- NCBI BLAST+
- HMMER 3.x
- RNAfold (ViennaRNA)
- Python 3.8+ (see requirements.txt)

