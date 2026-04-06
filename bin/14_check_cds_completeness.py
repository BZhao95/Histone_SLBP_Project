#!/usr/bin/env python3
"""
Description: Validates the completeness of Histone CDS sequences by checking 
             for canonical start (ATG) and stop (TAA, TAG, TGA) codons.
Classification: Complete, 5' partial, 3' partial, or Incomplete.
"""

import os
import argparse
from Bio import SeqIO

def get_args():
    parser = argparse.ArgumentParser(description="Check CDS integrity for start/stop codons.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file.")
    parser.add_argument("-o", "--output", required=True, help="Output TSV file for completeness status.")
    return parser.parse_args()

def check_start_codon(seq):
    """Returns True if sequence starts with ATG."""
    return seq.startswith("ATG")

def check_stop_codon(seq):
    """Returns True if sequence ends with TAA, TAG, or TGA."""
    return seq[-3:] in ["TAA", "TAG", "TGA"]

def main():
    args = get_args()
    
    # Status counters for the summary report
    stats = {"Complete": 0, "5'end partial": 0, "3'end partial": 0, "Incomplete": 0}

    print(f"[INFO] Analyzing sequence coompleteness in {args.input}...")

    with open(args.output, "w") as out_file:
        out_file.write("Protein_Accession\tStatus\n")
        
        for record in SeqIO.parse(args.input, "fasta"):
            # Normalize sequence: Uppercase and remove any whitespace/internal gaps
            sequence = str(record.seq).upper().replace(" ", "").replace("-", "")
            
            if len(sequence) < 3:
                continue

            has_start = check_start_codon(sequence)
            has_stop = check_stop_codon(sequence)
            
            # Classification
            if has_start and has_stop:
                status = "Complete"
            elif not has_start and has_stop:
                status = "5'end partial"
            elif has_start and not has_stop:
                status = "3'end partial"
            else:
                status = "Incomplete"

            # Parsing headers
            try:
                # Extracts the accession ID after the _cds_ tag
                record_id = record.id.split("|")[1].split("_cds_")[1]
            except IndexError:
                # Fallback if the header format differs
                record_id = record.id

            out_file.write(f"{record_id}\t{status}\n")
            stats[status] += 1

    # --- Final Summary ---

    print(f"[SUCCESS] Completeness results saved to {args.output}")

if __name__ == "__main__":
    main()
