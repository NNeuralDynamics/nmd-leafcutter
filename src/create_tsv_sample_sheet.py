#!/usr/bin/env python3
"""
Write one BAM list per cell type, no header, tab-separated:
    <bam_path>  <treatment>  <donor>

Parses names like:
    PS_STC1505_7_FB001VSmg1i_S7_L001.markdup.sorted.bam
    PS_STC1505_12_FB029TDMSOClean_S12_L001.markdup.sorted.bam
    cell type = FB, donor = 001V, treatment = Smg1i / DMSO

Usage:
    python make_bam_tsvs.py /path/to/bam_dir -o bam_lists/
    python make_bam_tsvs.py /path/to/bam_dir -o bam_lists/ --exclude DDALI DO
"""

import argparse
import os
import re
import sys
from collections import defaultdict

PATTERN = re.compile(
    r"(?:^|_)"                         # start of name, or after a prefix like PS_STC1505_7_
    r"(?P<cell_type>[A-Za-z]+\d?)"     # AT2, FB, MV, DD, DO, DDALI
    r"(?P<donor>\d{3}[A-Za-z])"        # 001V, 027U, 029T ...
    r"_?(?P<treatment>DMSO|Smg1i)"     # with or without an underscore before it
    r"(?:_?[Cc]lean)?"                 # _clean or Clean
    r"[._]"                            # followed by .markdup... or _S7_...
)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("bam_dir")
    p.add_argument("-o", "--outdir", default=".")
    p.add_argument("--exclude", nargs="+", default=["DDALI"],
                   help="Cell types to leave out (default: DDALI)")
    args = p.parse_args()

    bam_dir = os.path.abspath(args.bam_dir)
    groups, skipped = defaultdict(list), []

    for fname in sorted(os.listdir(bam_dir)):
        if not fname.endswith(".bam"):
            continue
        m = PATTERN.search(fname)
        if not m:
            skipped.append(fname)
            continue
        cell = m["cell_type"]
        if cell in args.exclude:
            continue
        path = os.path.join(bam_dir, fname)
        if not os.path.isfile(path + ".bai"):
            print(f"WARNING: no .bai index for {fname}", file=sys.stderr)
        groups[cell].append((path, m["treatment"], m["donor"]))

    if skipped:
        print("Did not match the naming pattern (skipped):", file=sys.stderr)
        for f in skipped:
            print(f"  {f}", file=sys.stderr)

    os.makedirs(args.outdir, exist_ok=True)
    for cell, rows in sorted(groups.items()):
        rows.sort(key=lambda r: (r[1] != "Smg1i", r[2]))  # Smg1i first, then by donor
        out = os.path.join(args.outdir, f"{cell}.txt")
        with open(out, "w") as fh:
            for r in rows:
                fh.write("\t".join(r) + "\n")
        donors = sorted({r[2] for r in rows})
        print(f"{out}: {len(rows)} BAMs, donors={','.join(donors)}")
        for d in donors:
            trts = {r[1] for r in rows if r[2] == d}
            if trts != {"DMSO", "Smg1i"}:
                print(f"  WARNING: {cell} {d} only has {', '.join(trts)}", file=sys.stderr)


if __name__ == "__main__":
    main()