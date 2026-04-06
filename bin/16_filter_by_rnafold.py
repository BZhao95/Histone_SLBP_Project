#!/usr/bin/env python3
"""
Description: Filters RNAfold results using the exact 3-line block parsing 
             and regex-based header extraction provided in the project methods.
Selection: Identifies the best (lowest MFE) canonical SL per accession.
"""

import re
import argparse
import os
from collections import defaultdict

def get_args():
    parser = argparse.ArgumentParser(description="Filter RNAfold output for best SL per gene.")
    parser.add_argument("-i", "--input", required=True, help="RNAfold output file (txt).")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA for deetected SLs.")
    return parser.parse_args()

def extract_and_find_best_sl(input_file, output_file):
    # Dictionary to group hits: { accession: [(position, mfe, header, sequence), ...] }
    sl_dict = defaultdict(list)
    
    # The dot-bracket structure (((((....)))))
    target_structure = r"\(\(\(\(\(\(\.\.\.\.\)\)\)\)\)"

    print(f"[INFO] Reading RNAfold results from {input_file}...")

    with open(input_file, "r") as infile:
        lines = infile.readlines()

        #3 lines per sequence
        for i in range(0, len(lines) - 2, 3):
            header = lines[i].strip()
            sequence = lines[i + 1].strip()
            rnafold_line = lines[i + 2].strip()

            # 1. Check length (26 bp) and verify the fold
            if len(sequence) == 26 and re.search(target_structure, rnafold_line):
                
                # 2. Extract accession and position
                # Matches format: >Accession/extended_START-END
                match = re.match(r">(\S+)/extended_(\d+)-\d+", header)
                if match:
                    accession = match.group(1)
                    position = int(match.group(2))

                    # 3. Extract MFE
                    #captures positive or negative floats
                    rnafold_match = re.search(r"\(([-+]?[0-9]*\.?[0-9]+)\)", rnafold_line)
                    rnafold_score = float(rnafold_match.group(1)) if rnafold_match else float("inf")

                    # Store for comparison
                    sl_dict[accession].append((position, rnafold_score, header, sequence))

    # 4. Filter for the "Best" SL per Accession
    print(f"[INFO] Filtering for best Minimum Free Energy (MFE)...")
    with open(output_file, "w") as outfile:
        count = 0
        for accession, entries in sl_dict.items():
            # Find entry with the lowest (most negative) RNAfold score
            best_entry = min(entries, key=lambda x: x[1])
            best_header, best_sequence = best_entry[2], best_entry[3]
            mfe_val = best_entry[1]
            
            outfile.write(f"{best_header} MFE={mfe_val}\n{best_sequence}\n")
            count += 1

    print(f"[SUCCESS] Filtered {count} unique best stem-loops into {output_file}")

if __name__ == "__main__":
    args = get_args() 
    extract_and_find_best_sl(args.input, args.output)
