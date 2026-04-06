#!/usr/bin/env python3
"""
Description: Generates a high-resolution presence/absence matrix for SLBP and co-factors.
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input meta.csv file")
    parser.add_argument("-o", "--output", default="presence_matrix.png", help="Output image path")
    return parser.parse_args()

def main():
    args = get_args()
    df = pd.read_csv(args.input, index_col=0)
    
    # Map values to symbols
    symbol_map = {0: '□', 1: '▨', 2: '■'}
    color_map = {0: '#8dc0ef', 1: '#073763', 2: '#054b89'}

    plt.figure(figsize=(15, 18))
    ax = plt.gca()

    x_spacing, y_spacing = 2.0, 1.0

    for y, (name, row) in enumerate(df.iterrows()):
        # Taxonomic Background Shading Logic
        bg_color = 'white'
        if 'Bilateria' in name: bg_color = '#F8F2FC'
        elif 'Fungi' in name: bg_color = '#FBF6E8'
        elif any(g in name for g in ['Chlorophyta', 'Streptophyta']): bg_color = '#F1FBE9'
        elif any(g in name for g in ['Alveolata', 'Stramenopiles']): bg_color = '#FBF0F0'

        ax.add_patch(plt.Rectangle((-0.5, y - 0.5), len(df.columns) * x_spacing, 1, 
                                   facecolor=bg_color, edgecolor='none', zorder=-1))

        for x, (col, val) in enumerate(row.items()):
            ax.text(x * x_spacing, y, symbol_map.get(val, '?'), ha='center', va='center',
                    fontsize=38, color=color_map.get(val, 'black'), fontfamily='DejaVu Sans')

    ax.set_yticks(np.arange(len(df)))
    ax.set_yticklabels(df.index, fontsize=12)
    ax.set_xticks(np.arange(0, len(df.columns) * 2, 2))
    ax.set_xticklabels(df.columns, fontsize=14, rotation=45)
    
    ax.invert_yaxis()
    plt.savefig(args.output, dpi=300, bbox_inches='tight')
    print(f"[SUCCESS] Matrix plot saved to {args.output}")

if __name__ == "__main__":
    main()
