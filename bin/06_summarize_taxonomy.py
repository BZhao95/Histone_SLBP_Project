#!/usr/bin/env python3
"""
Description: Summarizes ortholog hits across eukaryotic taxonomy using ETE3.
Features: 
- Filters for Eukaryota (TaxID 2759)
- Collapses major clades for better visualization
- Adds color-coded markers for different search sources (BLASTP, OrthoDB, etc.)
"""

import os
import argparse
import pandas as pd
from ete3 import NCBITaxa, TreeStyle, TextFace, CircleFace

def get_args():
    parser = argparse.ArgumentParser(description="Summarize taxonomic hits and render a tree.")
    parser.add_argument("-i", "--input_dir", required=True, help="Directory containing TaxID files.")
    parser.add_argument("-o", "--output_dir", required=True, help="Directory to save the results.")
    return parser.parse_args()

def is_eukaryota(taxid, ncbi):
    """Check if a taxid belongs to the Eukaryota (NCBI taxon: 2759)."""
    try:
        lineage = ncbi.get_lineage(int(taxid))
        return 2759 in lineage
    except:
        return False

def main():
    args = get_args()
    ncbi = NCBITaxa()
    os.makedirs(args.output_dir, exist_ok=True)

    # merge all taxids from different sources
    files = {
        'BLASTP': 'BLASTP_taxids.txt',
        'OrthoDB': 'OrthoDB_taxids.txt',
        'InterPro': 'InterPro_taxids.txt',
        'TBLASTN': 'tblastn/TBLASTN_taxids.txt'
    }

    # Load and combine data
    dfs = []
    for label, filename in files.items():
        path = os.path.join(args.input_dir, filename)
        if os.path.exists(path):
            tmp_df = pd.read_csv(path, sep='\t', header=None, names=['taxid'], dtype={'taxid': str})
            tmp_df['source'] = label
            dfs.append(tmp_df)
    
    if not dfs:
        print("[ERROR] No input TaxID files found in the specified directory.")
        return

    combined_df = pd.concat(dfs)
    combined_df['taxid'] = combined_df['taxid'].astype(str)

    # Filter for Eukaryota
    print("[INFO] Filtering for Eukaryotic lineages...")
    unique_ids = combined_df['taxid'].unique().tolist()
    filtered_ids = [tid for tid in unique_ids if is_eukaryota(tid, ncbi)]
    
    # Generate topology
    tree = ncbi.get_topology(filtered_ids, intermediate_nodes=True)
    name_map = ncbi.get_taxid_translator(filtered_ids)
    source_map = combined_df.groupby('taxid')['source'].apply(lambda x: set(x)).to_dict()

    # Mark the sources
    for node in tree.traverse():
        if node.is_leaf():
            taxid = node.name
            node.add_feature("original_taxid", str(taxid))
            node.name = name_map.get(int(taxid), taxid)
            
            sources = source_map.get(str(taxid), set())
            # Map sources to colors
            colors = {'BLASTP': 'red', 'OrthoDB': 'yellow', 'InterPro': 'blue', 'TBLASTN': 'green'}
            for idx, (src, color) in enumerate(colors.items()):
                if src in sources:
                    node.add_face(CircleFace(10, color), column=idx + 1, position='aligned')

    # Collapse
    collapse_taxa = ["Bilateria", "Cnidaria", "Placozoa", "Ctenophora", "Porifera", "Choanoflagellata",
                     "Filasterea", "Ichthyosporea", "Aphelida", "Fungi", "Amoebozoa", "Streptophyta",
                     "Chlorophyta", "Prasinodermophyta", "Stramenopiles", "Alveolata", "Rhizaria",
                     "Haptista", "Metamonada", "Discoba", "Rhodophyta", "Cryptophyceae"]

    for taxon in collapse_taxa:
        taxid_map = ncbi.get_name_translator([taxon])
        if taxon in taxid_map:
            target_id = str(taxid_map[taxon][0])
            nodes = tree.search_nodes(name=target_id)
            if nodes:
                n = nodes[0]
                num_leaves = len(n.get_leaves())
                n.children = [] # Collapse
                n.name = f"{taxon} ({num_leaves} species)"

    # Render
    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.scale = 50
    ts.title.add_face(TextFace("Phylogenetic Distribution Summary", fsize=20), column=0)
    
    output_png = os.path.join(args.output_dir, "taxonomy_summary.png")
    tree.render(output_png, w=2000, units='mm', tree_style=ts)
    tree.write(format=1, outfile=os.path.join(args.output_dir, "taxonomy_summary.nwk"))
    print(f"[SUCCESS] Results saved to {args.output_dir}")

if __name__ == "__main__":
    main()
