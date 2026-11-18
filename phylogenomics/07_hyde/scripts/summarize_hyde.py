#!/usr/bin/env python3
"""
Summarize the HyDe results.

Usage:
    python3 summarize_hyde.py hyde-out-filtered.txt
    python3 summarize_hyde.py hyde-out-filtered.txt hyde-ind.txt hyde_map.txt

1. Lists the significant triples (hyde-out-filtered.txt) sorted by Z-score.
2. Counts in how many significant triples each population is the hybrid.
3. (optional) For hyde-ind.txt (with the map file): counts, for each triple, how many
   individuals of the hybrid population are significant
   (p-value < 0.05 / number of tests, Bonferroni, and 0 < gamma < 1).
"""
import sys
import csv
from collections import Counter, defaultdict


def read_table(path):
    with open(path) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    for r in rows:
        r["Zscore"] = float(r["Zscore"])
        r["Pvalue"] = float(r["Pvalue"])
        r["Gamma"] = float(r["Gamma"])
    return rows


if len(sys.argv) < 2:
    sys.exit(__doc__)

# --- 1. significant triples ----------------------------------------------
triples = read_table(sys.argv[1])
triples.sort(key=lambda r: r["Zscore"], reverse=True)
print(f"## {len(triples)} significant triples (sorted by Z-score)\n")
print("P1\tHybrid\tP2\tZscore\tPvalue\tGamma")
for r in triples:
    print(f"{r['P1']}\t{r['Hybrid']}\t{r['P2']}\t{r['Zscore']:.2f}\t{r['Pvalue']:.2e}\t{r['Gamma']:.3f}")

# --- 2. how often each population is the hybrid --------------------------
print("\n## Number of significant triples per hybrid population\n")
for pop, n in Counter(r["Hybrid"] for r in triples).most_common():
    print(f"{pop}\t{n}")

# --- 3. individuals (optional) -------------------------------------------
if len(sys.argv) > 3:
    ind = read_table(sys.argv[2])
    pop_of = {}
    with open(sys.argv[3]) as f:
        for line in f:
            if line.strip():
                sample, pop = line.split()[:2]
                pop_of[sample] = pop
    cutoff = 0.05 / len(ind)
    print(f"\n## Individuals (Bonferroni: p < 0.05 / {len(ind)} = {cutoff:.2e})\n")
    print("P1\tHybrid_pop\tP2\tsignificant/total\tmean_gamma_significant")
    groups = defaultdict(list)
    for r in ind:
        groups[(r["P1"], pop_of.get(r["Hybrid"], "?"), r["P2"])].append(r)
    for (p1, hyb, p2), rows in groups.items():
        sig = [r for r in rows if r["Pvalue"] < cutoff and 0 < r["Gamma"] < 1]
        mean_g = sum(r["Gamma"] for r in sig) / len(sig) if sig else float("nan")
        print(f"{p1}\t{hyb}\t{p2}\t{len(sig)}/{len(rows)}\t{mean_g:.3f}")
elif len(sys.argv) == 3:
    print("\nTo summarize hyde-ind.txt, give the map file too (see usage).")
