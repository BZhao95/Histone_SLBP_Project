#!/bin/bash

# ==============================================================================
# Script: 13_get_flanking_sequences.sh
# Description: Retrieves 300bp downstream of histone stop codons.
#              Handles strand orientation and genomic coordinate math.
# ==============================================================================

# --- Settings ---
PADDING=300
INPUT_TSV="results/histones/meta_info.tsv" 
OUT_DIR="data/flanking_sequences"
LOG_DIR="logs"

mkdir -p "$OUT_DIR" "$LOG_DIR"
UPSTREAM_ISSUES="$LOG_DIR/upstream_issues.txt"
DOWNSTREAM_ISSUES="$LOG_DIR/downstream_issues.txt"


# --- Logic Loop ---
# Skips the header and reads metadata
tail -n +2 "$INPUT_TSV" | while IFS=$'\t' read -r prot_acc nuc_acc cdd cdd_name location completeness; do
    
    # Skip incomplete sequences per Methods (Marzluff et al. 2008 requirement)
    if [[ "$completeness" != "Complete" ]]; then
        continue
    fi

    echo "Processing $prot_acc..."

    # Normalize location string (remove 'complement' and 'join' wrappers)
    clean_loc=$(echo "$location" | sed -E 's/complement\(([^)]+)\)/\1/g; s/join\(([^)]+)\)/\1/g')
    
    # Extract genomic boundaries
    start=$(echo "$clean_loc" | grep -oE '[0-9]+' | head -n 1)
    end=$(echo "$clean_loc" | grep -oE '[0-9]+' | tail -n 1)

    # --- Coordinate Calculation ---
    if [[ $location == complement* ]]; then
        # COMPLEMENT STRAND: 3' end is at the lower coordinate
        stop_pos=$start
        ext_start=$((stop_pos + 2))        # Includes end of CDS for SL context
        ext_end=$((stop_pos - PADDING))    # Moves "down" the genome
        
        if ((ext_end < 1)); then
            echo "$prot_acc" >> "$DOWNSTREAM_ISSUES"
            continue
        fi
    else
        # FORWARD STRAND: 3' end is at the higher coordinate
        stop_pos=$end
        ext_start=$((stop_pos - 2))        # Includes end of CDS
        ext_end=$((stop_pos + PADDING))    # Moves "up" the genome
        
        if ((ext_start < 1)); then
            echo "$prot_acc" >> "$UPSTREAM_ISSUES"
            ext_start=1
        fi
    fi

    # --- Fetching ---
    TEMP_FASTA="$OUT_DIR/${prot_acc}_temp.fasta"
    FINAL_FASTA="$OUT_DIR/${prot_acc}_downstream.fasta"

    if ! efetch -db nucleotide -id "$nuc_acc" -seq_start "$ext_start" -seq_stop "$ext_end" -format fasta > "$TEMP_FASTA"; then
        echo "[ERROR] Failed to fetch $prot_acc"
        continue
    fi

    # Clean FASTA header to match Protein ID
    sed -i.bak "1s/.*/>${prot_acc}/" "$TEMP_FASTA" && rm "${TEMP_FASTA}.bak"
    mv "$TEMP_FASTA" "$FINAL_FASTA"

    sleep 1 # Be polite to NCBI servers
done

echo "[FINISH] Downstream flanking retrieval complete."
