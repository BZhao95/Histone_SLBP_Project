#!/usr/bin/env python3
"""
Description: Scans histone 3' flanking sequences for canonical and 
             non-canonical Polyadenylation Signals (PAS).
Logic: Uses a sliding regex search to identify all occurrences of 15 
       PAS variants within the 300bp downstream region.
Output: TSV file containing protein ID, the specific motif found, and its position.
"""

import argparse
import re
import os
from Bio import SeqIO

def get_args():
    parser = argparse.ArgumentParser(description="Scan for PAS motifs in histone flanking regions.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA (300bp flanking sequences).")
    parser.add_argument("-o", "--output", required=True, help="Output TSV file path.")
    return parser.parse_args()

# Canonical and frequent PAS variants 
PAS_VARIANTS = [
    "AATAAA", "ATTAAA", "TATAAA", "AAGAAA", "AGTAAA",
    "AATACA", "AATATA", "AATAAG", "CATAAA", "GATAAA",
    "AAATAA", "AAATTA", "TAATAA", "AAAATA", "AATAAT"
]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fasta", help="FASTA file with flanking sequences")
    ap.add_argument("--out", default="pas_hits.tsv", help="Output TSV")
    args = ap.parse_args()

    with open(args.out, "w") as out:
        out.write("Protein_Accession\tPAS_motif\tposition\n")
        for record in SeqIO.parse(args.fasta, "fasta"):
            sid = record.id
            seq = str(record.seq).upper()
            for motif in PAS_VARIANTS:
                for m in re.finditer(motif, seq):
                    pos = m.start() + 1  # 1-based
                    out.write(f"{sid}\t{motif}\t{pos}\n")
    print(f"[+] PAS search finished. Results saved to {args.out}")

if __name__ == "__main__":
    main()
