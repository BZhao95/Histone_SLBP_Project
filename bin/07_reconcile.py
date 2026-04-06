#!/usr/bin/env python3
"""
Description: Reconciles an SLBP gene tree with a eukaryotic species tree from NCBI taxonomy.
Algorithm: Least Common Ancestor (LCA) mapping.
Outputs: Annotated Newick tree and frequency summary of events.
"""

import argparse
from ete3 import Tree, NCBITaxa

def get_args():
    parser = argparse.ArgumentParser(description="Reconcile Gene Tree with Species Tree.") 
    parser.add_argument("-g", "--gene_tree", required=True, help="Input gene tree (Newick)") #FastTree gene tree
    parser.add_argument("-s", "--species_tree", required=True, help="Input species tree (Newick)") #directly from NCBI Taxonomy
    parser.add_argument("-m", "--map", required=True, help="Gene-to-Species mapping (TSV: gene_id\tSpecies_Name)") # gene_id\tSpecies_name_with_underscores
    parser.add_argument("-o", "--out", required=True, help="Output prefix")
    return parser.parse_args()

def main():
    args = get_args()
    ncbi = NCBITaxa()

    print("[INFO] Loading trees...")
    gt = Tree(args.gene_tree, format=1)
    st = Tree(args.species_tree, format=1)
    species_tree_taxids = {leaf.name for leaf in st.iter_leaves()}

    # Map gene IDs to TaxIDs
    gene2taxid = {}
    with open(args.map) as f:
        for line in f:
            if line.startswith("#") or not line.strip(): continue
            gid, sp = line.strip().split("\t")[:2]
            sp_name = sp.replace("_", " ")
            try:
                trans = ncbi.get_name_translator([sp_name])
                if sp_name in trans:
                    gene2taxid[gid] = str(trans[sp_name][0])
            except: pass

    # LCA Mapping
    for node in gt.traverse("postorder"):
        if node.is_leaf():
            tid = gene2taxid.get(node.name)
            node.add_features(species=tid, species_set={tid} if tid else set())
        else:
            sp_set = set()
            for ch in node.children:
                sp_set |= getattr(ch, "species_set", set())
            node.species_set = sp_set
            
            if sp_set:
                try:
                    lca = st.get_common_ancestor(list(sp_set))
                    node.add_features(sp_lca=lca.name)
                except: node.add_features(sp_lca=None)

    # Classification
    dups, specs = [], []
    for node in gt.traverse("postorder"):
        if node.is_leaf(): 
            node.add_features(event="leaf")
            continue
        
        child_lcas = {getattr(ch, "sp_lca", None) for ch in node.children}
        if node.sp_lca in child_lcas:
            node.add_features(event="duplication")
            dups.append(node)
        else:
            node.add_features(event="speciation")
            specs.append(node)

    # Save results
    gt.write(format=1, outfile=f"{args.out}_reconciled.nwk")
    with open(f"{args.out}_summary.tsv", "w") as out:
        out.write(f"duplication\t{len(dups)}\nspeciation\t{len(specs)}\n")
    
    print(f"[OUTPUT] Reconciled tree saved. Duplications: {len(dups)}, Speciations: {len(specs)}")

if __name__ == "__main__":
    main()
