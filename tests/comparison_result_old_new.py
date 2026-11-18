#!/usr/bin/env python3
"""
Compare two LeafCutter DS runs of the same data, e.g. treatment only
(old) vs treatment + donor covariate (new).

Expects in each run directory:  <run_dir>/<CT>/<CT>_cluster_significance.txt
                                <run_dir>/<CT>/<CT>_effect_sizes.txt

Cluster IDs (clu_N) are renumbered on every run, so clusters are matched
by the intron coordinates they contain, not by ID.

Usage:
    python compare_leafcutter_runs.py \
        --old_dir /data/shared/sthakur/nmd-FB_MV/leafcutter_output_final_RF \
        --new_dir /scratch/thakur.san/nmd-leafcutter/leafcutter_donor_results_v1 \
        --cell_types AT2 DD DO FB MV \
        --out comparison_out
"""

import argparse
import os
from collections import Counter

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load_run(run_dir, ct):
    base = os.path.join(run_dir, ct, ct)
    sig = pd.read_csv(f"{base}_cluster_significance.txt", sep="\t")
    eff = pd.read_csv(f"{base}_effect_sizes.txt", sep="\t")

    parts = eff["intron"].str.split(":")
    strand = parts.str[3].str[-1]                       # clu_123_+  ->  +
    eff["coord"] = parts.str[0] + ":" + parts.str[1] + ":" + parts.str[2] + ":" + strand
    eff["cluster"] = parts.str[0] + ":" + parts.str[3]
    eff = eff.drop_duplicates("coord")  

    # largest |deltapsi| per cluster, keeping its sign
    idx = eff.groupby("cluster")["deltapsi"].apply(lambda s: s.abs().idxmax())
    max_dpsi = eff.loc[idx.values, ["cluster", "deltapsi"]].set_index("cluster")["deltapsi"]
    sig = sig.set_index("cluster")
    sig["max_dpsi"] = max_dpsi
    return sig, eff


def match_clusters(eff_old, eff_new):
    """old cluster -> new cluster sharing the most intron coordinates."""
    coord_to_new = dict(zip(eff_new["coord"], eff_new["cluster"]))
    mapping = {}
    for clu, coords in eff_old.groupby("cluster")["coord"]:
        hits = Counter(coord_to_new[c] for c in coords if c in coord_to_new)
        if hits:
            mapping[clu] = hits.most_common(1)[0][0]
    return mapping


def is_sig(df, sfx, fdr, min_dpsi):
    return ((df[f"status{sfx}"] == "Success") & (df[f"p.adjust{sfx}"] < fdr)
            & (df[f"max_dpsi{sfx}"].abs() >= min_dpsi))


def compare(ct, old_dir, new_dir, out_dir, fdr, min_dpsi):
    sig_o, eff_o = load_run(old_dir, ct)
    sig_n, eff_n = load_run(new_dir, ct)
    mapping = match_clusters(eff_o, eff_n)

    m = pd.DataFrame({"cluster_old": list(mapping), "cluster_new": list(mapping.values())})
    m = (m.join(sig_o.add_suffix("_old"), on="cluster_old")
          .join(sig_n.add_suffix("_new"), on="cluster_new"))
    m = m[(m["status_old"] == "Success") & (m["status_new"] == "Success")].copy()

    so = is_sig(m, "_old", fdr, min_dpsi)
    sn = is_sig(m, "_new", fdr, min_dpsi)
    m["category"] = np.select([so & sn, ~so & sn, so & ~sn], ["both", "new only", "old only"], "neither")

    # intron-level deltapsi agreement
    d = pd.concat([eff_o.set_index("coord")["deltapsi"], eff_n.set_index("coord")["deltapsi"]],
                  axis=1, keys=["old", "new"]).dropna()

    lp_o = -np.log10(m["p_old"].clip(lower=1e-300))
    lp_n = -np.log10(m["p_new"].clip(lower=1e-300))
    either = m[m["category"] != "neither"]
    same_dir = (np.sign(either["max_dpsi_old"]) == np.sign(either["max_dpsi_new"])).mean()

    summary = {
        "cell_type": ct,
        "tested_old": int((sig_o["status"] == "Success").sum()),
        "tested_new": int((sig_n["status"] == "Success").sum()),
        "matched_tested_both": len(m),
        "sig_old": int(is_sig(sig_o, "", fdr, min_dpsi).sum()),
        "sig_new": int(is_sig(sig_n, "", fdr, min_dpsi).sum()),
        "sig_both": int((m["category"] == "both").sum()),
        "gained_with_donor": int((m["category"] == "new only").sum()),
        "lost_with_donor": int((m["category"] == "old only").sum()),
        "spearman_log10p": round(pd.Series(lp_o.values).corr(pd.Series(lp_n.values), method="spearman"), 3),
        "pearson_deltapsi": round(d["old"].corr(d["new"]), 3),
        "same_direction_sig": round(same_dir, 3) if len(either) else np.nan,
    }

    cols = ["cluster_old", "cluster_new", "genes_new", "p.adjust_old", "p.adjust_new",
            "max_dpsi_old", "max_dpsi_new"]
    for cat, name in [("new only", "gained"), ("old only", "lost")]:
        (m.loc[m["category"] == cat, cols].sort_values("p.adjust_new")
          .to_csv(os.path.join(out_dir, f"{ct}_{name}_with_donor.tsv"), sep="\t", index=False))

    # plots
    colors = {"both": "#2b6cb0", "new only": "#d69e2e", "old only": "#c53030", "neither": "#bbbbbb"}
    fig, ax = plt.subplots(1, 2, figsize=(11, 5))
    for cat in ["neither", "both", "old only", "new only"]:
        sel = m["category"] == cat
        ax[0].scatter(lp_o[sel], lp_n[sel], s=4, alpha=0.5, c=colors[cat],
                      label=f"{cat} ({sel.sum()})", rasterized=True)
    lim = max(lp_o.max(), lp_n.max())
    ax[0].plot([0, lim], [0, lim], "k--", lw=0.8)
    ax[0].set_xlabel("-log10 p, treatment only")
    ax[0].set_ylabel("-log10 p, treatment + donor")
    ax[0].set_title(f"{ct}: cluster p-values")
    ax[0].legend(markerscale=3, fontsize=8, title=f"FDR < {fdr}", title_fontsize=8)

    ax[1].scatter(d["old"], d["new"], s=3, alpha=0.3, c="#4a5568", rasterized=True)
    ax[1].plot([-1, 1], [-1, 1], "k--", lw=0.8)
    ax[1].set_xlim(-1, 1); ax[1].set_ylim(-1, 1)
    ax[1].set_xlabel("deltapsi, treatment only")
    ax[1].set_ylabel("deltapsi, treatment + donor")
    ax[1].set_title(f"{ct}: intron deltapsi (r = {summary['pearson_deltapsi']})")
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, f"{ct}_comparison.png"), dpi=200)
    plt.close(fig)
    return summary


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--old_dir", required=True, help="Run without donor covariate")
    p.add_argument("--new_dir", required=True, help="Run with donor covariate")
    p.add_argument("--cell_types", nargs="+", default=["AT2", "DD", "DO", "FB", "MV"])
    p.add_argument("--out", default="comparison_out")
    p.add_argument("--fdr", type=float, default=0.05)
    p.add_argument("--min_dpsi", type=float, default=0.0,
                   help="Also require |deltapsi| >= this for a hit (e.g. 0.1)")
    args = p.parse_args()

    os.makedirs(args.out, exist_ok=True)
    rows = []
    for ct in args.cell_types:
        try:
            rows.append(compare(ct, args.old_dir, args.new_dir, args.out, args.fdr, args.min_dpsi))
            print(f"{ct}: done")
        except FileNotFoundError as e:
            print(f"{ct}: skipped, missing file: {e.filename}")

    if rows:
        summary = pd.DataFrame(rows)
        summary.to_csv(os.path.join(args.out, "comparison_summary.tsv"), sep="\t", index=False)
        print()
        print(summary.to_string(index=False))


if __name__ == "__main__":
    main()