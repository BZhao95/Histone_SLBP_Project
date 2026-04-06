#!/usr/bin/env python3
import argparse
from ete3 import NCBITaxa

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--slbp", required=True)
    parser.add_argument("--lsm10", required=True)
    parser.add_argument("--lsm11", required=True)
    parser.add_argument("-q", "--query", required=True, help="Clade to search within (e.g. 'Metazoa')")
    args = parser.parse_args()

    ncbi = NCBITaxa()
    
    def load(p): return set(line.strip() for line in open(p))
    
    
    common = load(args.slbp) & load(args.lsm10) & load(args.lsm11)
    
    # Find ancestors
    query_id = ncbi.get_name_translator([args.query]).get(args.query, [None])[0]
    if query_id:
        descendants = ncbi.get_descendant_taxa(query_id, intermediate_nodes=True)
        matching = [tid for tid in descendants if str(tid) in common]
        names = ncbi.get_taxid_translator(matching)
        print(f"\nSpecies in {args.query} with complete SLBP/LSM10/LSM11 set:")
        for tid in matching:
            print(f"- {names[tid]} (TaxID: {tid})")

if __name__ == "__main__":
    main()
