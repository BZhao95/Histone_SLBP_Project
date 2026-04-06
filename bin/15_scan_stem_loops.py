#!/usr/bin/env python3
"""
Description: Scans histone 3' flanking regions for a 6-4-6 stem-loop structure.
Criteria: 6bp stem, 4nt loop, 6bp stem. Allows at most one mismatch.
Output: Extracted 16nt cores and 26nt extended sequences (for RNAfold).
"""

import os
import argparse
from Bio import SeqIO

def get_args():
    parser = argparse.ArgumentParser(description="Scan for histone stem-loops.")
    parser.add_argument("-i", "--input", required=True, help="Input flanking FASTA.")
    parser.add_argument("-o", "--output_dir", required=True, help="Directory for outputs.")
    return parser.parse_args()

def is_valid_sl(candidate):
    """Checks for a 6-4-6 structure with at most one mismatch."""
    if len(candidate) != 16: return False
    
    stem1 = candidate[:6].upper().replace('U', 'T')
    stem2 = candidate[10:].upper().replace('U', 'T')
    
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    mismatches = 0
    
    # Check Watson-Crick pairing (allow at most 1 pair of mismatch)
    for i in range(6):
        if complement.get(stem1[i]) != stem2[5-i]:
            mismatches += 1
            if mismatches > 1:
                return False
    return True

def main():
    args = get_args()
    os.makedirs(args.output_dir, exist_ok=True)
    
    out_detected = os.path.join(args.output_dir, "detected_cores.fasta")
    out_extended = os.path.join(args.output_dir, "extended_sequences.fasta")
    out_none = os.path.join(args.output_dir, "no_sl_found.fasta")

    counts = {"total": 0, "with_sl": 0}

    with open(out_detected, "w") as f_core, \
         open(out_extended, "w") as f_ext, \
         open(out_none, "w") as f_none:

        for record in SeqIO.parse(args.input, "fasta"):
            counts["total"] += 1
            seq = str(record.seq).upper()
            found_in_this_record = False

            for i in range(len(seq) - 16 + 1):
                window = seq[i:i+16]
                
                if is_valid_sl(window):
                    found_in_this_record = True
                    # Core 16nt
                    f_core.write(f">{record.id}_pos_{i+1}\n{window}\n")
                    
                    # Extended 26nt (add 5nt padding each side)
                    start = max(0, i - 5)
                    end = min(len(seq), i + 21)
                    ext_seq = seq[start:end]
                    if len(ext_seq) == 26:
                        f_ext.write(f">{record.id}_ext_{start+1}\n{ext_seq}\n")

            if found_in_this_record:
                counts["with_sl"] += 1
            else:
                f_none.write(f">{record.id}\n{seq}\n")

    print(f"[INFO] Processed {counts['total']} sequences.")
    print(f"[INFO] Found SL candidates in {counts['with_sl']} sequences.")

if __name__ == "__main__":
    main()
