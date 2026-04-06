#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.stats import gaussian_kde
from matplotlib.patches import Patch
import argparse
import sys

def main():
    parser = argparse.ArgumentParser(description="Compare processing states across basal fungal clades.")
    parser.add_argument("-i", "--input", help="Path to meta_fungi.tsv", 
                        default="data/histone_fungi/fungi_meta.tsv")
    args = parser.parse_args()

    tsv_path = args.input

    # file & column names
    COL_ACC   = "Protein_Accession"  
    COL_SPEC  = "Species"            
    COL_CLADE = "Clade"              
    COL_SLBP  = "SLBP"               
    COL_SL    = "SL"                 
    COL_PAS   = "PAS"                

    CLADES = ["Cryptomycota", "Chytridiomycota", "Blastocladiomycota", "Zoopagomycota", "Mucoromycota", "Microsporidia"]

    STATE_ORDER = ["SLBP–SL", "SLBP–SL–PAS", "SLBP–PAS", "PAS-only"]
    STATE_COLORS = {
        "SLBP–SL":      "#2ca02c",  
        "SLBP–SL–PAS":  "#1f77b4",  
        "SLBP–PAS":     "#ff7f0e",  
        "PAS-only":     "#d62728",  
    }

    # Load & basic cleaning
    df = pd.read_csv(tsv_path, sep="\t", dtype=str).fillna("NA")

    def is_yes(x):
        return str(x).strip().upper() == "Y"

    def has_pas(x):
        return (isinstance(x, str) and x.upper() != "NA" and len(x.strip()) > 0)

    df["_SLBP"] = df[COL_SLBP].apply(is_yes)
    df["_SL"]   = df[COL_SL].apply(is_yes)
    df["_PAS"]  = df[COL_PAS].apply(has_pas)

    # Classify locus state
    def classify_row(row):
        slbp, sl, pas = row["_SLBP"], row["_SL"], row["_PAS"]
        if slbp and sl and not pas:
            return "SLBP–SL"
        if slbp and sl and pas:
            return "SLBP–SL–PAS"
        if slbp and (not sl) and pas:
            return "SLBP–PAS"
        if (not slbp) and (not sl) and pas:
            return "PAS-only"
        return "Other"

    df["State"] = df.apply(classify_row, axis=1)
    df_focus = df[df["State"].isin(STATE_ORDER)].copy()

    # Metrics
    total_records = len(df)
    excluded_records = (df["State"] == "Other").sum()
    included_records = len(df_focus)

    print(f"Total records: {total_records}")
    print(f"Included in state analysis: {included_records}")
    print(f"Excluded (Other): {excluded_records}")
    print(f"Excluded percentage: {excluded_records/total_records*100:.2f}%")

    # Per-species composition
    def species_fraction(group):
        n = len(group)
        fr = group["State"].value_counts(normalize=True)
        d = {s: fr.get(s, 0.0) for s in STATE_ORDER}
        d["N_loci"] = n
        return pd.Series(d)

    species_comp = (
        df_focus
        .groupby([COL_SPEC, COL_CLADE])
        .apply(species_fraction)
        .reset_index()
        .rename(columns={COL_SPEC: "Species", COL_CLADE: "Clade"})
    )

    # Per-clade mean composition
    clade_mean = (
        species_comp
        .groupby("Clade")[STATE_ORDER]
        .mean()
        .reindex(CLADES)
        .fillna(0.0)
    )

    # Dominant state per species
    def dominant_state_row(row):
        vals = [row[s] for s in STATE_ORDER]
        if sum(vals) == 0:
            return "Other"
        mx = max(vals)
        ties = [STATE_ORDER[i] for i, v in enumerate(vals) if v == mx]
        return "SLBP–SL–PAS" if "SLBP–SL–PAS" in ties else ties[0]

    species_comp["Dominant_State"] = species_comp.apply(dominant_state_row, axis=1)
    dom_counts = (
        species_comp[species_comp["Dominant_State"].isin(STATE_ORDER)]
        .groupby(["Clade", "Dominant_State"])
        .size()
        .unstack(fill_value=0)
        .reindex(index=CLADES, columns=STATE_ORDER)
        .fillna(0)
    )

    # Log-odds contrast
    rng = np.random.default_rng(42)

    def safe_logit(p):
        eps = 1e-6
        p = np.clip(p, eps, 1-eps)
        return np.log(p / (1-p))

    species_comp["_p_SL"]      = (species_comp["SLBP–SL"] + species_comp["SLBP–SL–PAS"]).astype(float)
    species_comp["_p_PASonly"] = species_comp["PAS-only"].astype(float)

    def bootstrap_mean_lor(vals_sl, vals_pas, B=3000):
        lor_per_species = safe_logit(vals_sl) - safe_logit(vals_pas)
        if len(lor_per_species) == 0:
            return np.nan, (np.nan, np.nan)
        idx = rng.integers(0, len(lor_per_species), size=(B, len(lor_per_species)))
        boot = lor_per_species[idx].mean(axis=1)
        mean = float(np.mean(boot))
        lo, hi = np.percentile(boot, [2.5, 97.5])
        return mean, (float(lo), float(hi))

    rows = []
    for cl in CLADES:
        sub = species_comp[species_comp["Clade"] == cl]
        mean, (lo, hi) = bootstrap_mean_lor(sub["_p_SL"].values, sub["_p_PASonly"].values, B=3000)
        rows.append({"Clade": cl, "mean_LOR": mean, "CI_low": lo, "CI_high": hi})
    lor_df = pd.DataFrame(rows).set_index("Clade").reindex(CLADES)

    # Save summaries
    out_species_tsv = tsv_path.replace(".tsv", "_species_state_fractions_3panels.tsv")
    out_clade_tsv   = tsv_path.replace(".tsv", "_clade_mean_state_fractions_3panels.tsv")
    species_comp.to_csv(out_species_tsv, sep="\t", index=False)
    clade_mean.to_csv(out_clade_tsv, sep="\t")

    # Figure Layout
    fig = plt.figure(figsize=(14, 10))
    gs = gridspec.GridSpec(2, 2, height_ratios=[1.2, 1])

    axA = fig.add_subplot(gs[0, :])  
    axB = fig.add_subplot(gs[1, 0])  
    axC = fig.add_subplot(gs[1, 1])  

    # Panel A: Ridgeline
    clade_y_bases = {cl: i*1.5 for i, cl in enumerate(reversed(CLADES))}
    x_eval = np.linspace(0, 1, 200)

    for cl in CLADES:
        sub = species_comp[species_comp["Clade"] == cl]
        y_base = clade_y_bases[cl]
        axA.plot([0, 1], [y_base, y_base], color='black', linewidth=0.8, alpha=0.5)
        
        for s in STATE_ORDER:
            vals = sub[s].values
            if len(vals) > 0 and vals.max() > 0:
                noise = np.random.normal(0, 1e-4, size=len(vals))
                vals_jittered = np.clip(vals + noise, 0, 1)
                try:
                    kde = gaussian_kde(vals_jittered, bw_method=0.1)
                    y_dens = kde(x_eval)
                    y_dens = (y_dens / max(y_dens.max(), 1)) * 1.2 
                    axA.fill_between(x_eval, y_base, y_base + y_dens, color=STATE_COLORS[s], alpha=0.6, zorder=y_base)
                    axA.plot(x_eval, y_base + y_dens, color=STATE_COLORS[s], linewidth=1.2, zorder=y_base)
                except:
                    axA.bar(vals.mean(), 1.2, bottom=y_base, width=0.01, color=STATE_COLORS[s], alpha=0.6, zorder=y_base)

    axA.set_yticks([clade_y_bases[cl] for cl in CLADES])
    axA.set_yticklabels(CLADES, fontsize=11)
    axA.set_xlim(0, 1)
    axA.set_xlabel("Species fraction within Clade", fontsize=11)
    axA.set_title("A", loc="left", fontsize=20)
    axA.spines['top'].set_visible(False)
    axA.spines['right'].set_visible(False)
    axA.spines['left'].set_visible(False)
    axA.tick_params(axis='y', length=0)
    legend_elements = [Patch(facecolor=STATE_COLORS[s], label=s, alpha=0.7) for s in STATE_ORDER]
    axA.legend(handles=legend_elements, frameon=False, fontsize=10, loc="upper right", bbox_to_anchor=(1.0, 1.05))

    # Panel B
    width = 0.14
    xB_axis = np.arange(len(dom_counts))
    for i, s in enumerate(STATE_ORDER):
        axB.bar(xB_axis + i*width, dom_counts[s].values, width=width, label=s, color=STATE_COLORS[s])
    axB.set_xticks(xB_axis + width*(len(STATE_ORDER)-1)/2)
    axB.set_xticklabels(dom_counts.index, rotation=35, ha="center")
    axB.set_ylabel("Species count")
    axB.set_title("B", loc="left", fontsize=20)

    # Panel C
    xC_axis = np.arange(len(lor_df))
    y = lor_df["mean_LOR"].values
    yerr = np.vstack([y - lor_df["CI_low"].values, lor_df["CI_high"].values - y])
    axC.errorbar(xC_axis, y, yerr=yerr, fmt='o', capsize=4, color='k')
    axC.axhline(0, linestyle="--", linewidth=1, color="grey")
    axC.set_xticks(xC_axis)
    axC.set_xticklabels(lor_df.index, rotation=35, ha="center")
    axC.set_ylabel("Log-odds: (SL-containing) vs (PAS-only)")
    axC.set_title("C", loc="left", fontsize=20)

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.1, hspace=0.35) 
    out_fig = tsv_path.replace(".tsv", "_3_panels_ridgeline.png")
    plt.savefig(out_fig, dpi=300, bbox_inches="tight")
    print(f"Saved 3-panel figure: {out_fig}")

if __name__ == "__main__":
    main()
