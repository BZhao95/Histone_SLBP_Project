#!/usr/bin/env python3
"""
Description: Scans the N-terminal region of SLBP sequences for the SLIP1-interacting motif.
Reference Motif: ARCKDWGSAVEEDEQL (from Xenopus laevis SLBP1)
"""

import os
import csv
import argparse
from Bio import SeqIO
from Bio import pairwise2

def get_args():
    parser = argparse.ArgumentParser(description="Scan SLBP sequences for SLIP1-binding motifs.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file of SLBP sequences.")
    parser.add_argument("-o", "--output_dir", required=True, help="Directory to save TSV and motif FASTA.")
    parser.add_argument("-s", "--score_min", type=float, default=20.0, help="Minimum alignment score (default: 20).")
    parser.add_argument("-id", "--identity_min", type=float, default=40.0, help="Minimum identity percentage (default: 40).")
    return parser.parse_args()

def main():
    args = get_args()
    
    # Reference motif from Xenopus laevis SLBP1 from the paper: SLIP1, a Factor Required for...
    REF_MOTIF = "ARCKDWGSAVEEDEQL"
    
    os.makedirs(args.output_dir, exist_ok=True)
    tsv_path = os.path.join(args.output_dir, "slip1_motif_hits.tsv")
    fasta_path = os.path.join(args.output_dir, "slip1_motifs_extracted.fasta")

    # Open output files
    with open(tsv_path, "w", newline="") as out_tsv, open(fasta_path, "w") as out_fasta:
        writer = csv.writer(out_tsv, delimiter="\t")
        writer.writerow([
            "Protein_ID", "Score", "Percent_Identity", "Hit_Start", "Hit_End", "Matched_Region"
        ])

        print(f"[INFO] Scanning {args.input} for SLIP1-interacting motifs...")
        hit_count = 0

        for record in SeqIO.parse(args.input, "fasta"):
            # Per methods: Extract the N-terminal 150 amino acids
            protein_id = record.id
            protein_seq = str(record.seq[:150])

            # Local alignment using standard parameters (Match: 2, Mismatch: -1, GapOpen: -0.5, GapExt: -0.1)
            alignments = pairwise2.align.localms(
                REF_MOTIF, protein_seq, 2, -1, -0.5, -0.1
            )

            if not alignments:
                continue

            # Take the best alignment
            best_align = alignments[0]
            score = best_align[2]
            aligned_motif = best_align[0]
            aligned_protein = best_align[1]
            # Coordinates are 0-indexed
            align_start = best_align[3]
            align_end = best_align[4]

            # Calculate Percent Identity
            matches = sum(m == p for m, p in zip(aligned_motif, aligned_protein) if m != '-' and p != '-')
            length = sum(1 for m, p in zip(aligned_motif, aligned_protein) if m != '-' and p != '-')
            identity = (matches / length) * 100 if length > 0 else 0

            # Filter based on thresholds defined in Methods
            if score >= args.score_min and identity >= args.identity_min:
                # Extract the specific motif sequence from the original protein
                matched_clean = protein_seq[align_start:align_end]
                
                # Save to TSV
                writer.writerow([
                    protein_id, f"{score:.1f}", f"{identity:.1f}",
                    align_start, align_end, matched_clean
                ])
                
                # Save to FASTA for downstream motif logo analysis (e.g., WebLogo)
                out_fasta.write(f">{protein_id}_motif\n{matched_clean}\n")
                hit_count += 1

    print(f"[SUCCESS] Processed sequences. Found {hit_count} motif-positive SLBPs.")
    print(f"[INFO] Results saved to {args.output_dir}")

if __name__ == "__main__":
    main()
