#!/usr/bin/env python3
"""
Summarize Dsuite Dtrios output.

Works with sets named by population (e.g. PopA) or by individual, named
<Population>_<ID> (e.g. PopA_1023): the trailing _<number> is removed to get
the population.

For each Dtrios file (_tree, _BBAA, _Dmin):
  - corrects the p-values for multiple testing (Benjamini-Hochberg and Bonferroni);
  - counts the significant trios (BH q < 0.05) by the populations of P2 and P3
    (a positive D means excess allele sharing between P2 and P3);
  - lists the most significant trios.

Usage:
    python3 summarize_dsuite.py <SETS prefix> [alpha]

Example:
    python3 summarize_dsuite.py my_SETS
    → reads my_SETS_tree.txt, _BBAA.txt and _Dmin.txt
    → writes my_SETS_<kind>_BH.txt (all trios with q-values)
"""
import csv
import os
import re
import sys
from collections import Counter

if len(sys.argv) < 2:
    sys.exit(__doc__)
prefix = sys.argv[1]
alpha = float(sys.argv[2]) if len(sys.argv) > 2 else 0.05


def pop(name):
    return re.sub(r'_\d+$', '', name)


def bh(pvalues):
    """Benjamini-Hochberg q-values."""
    n = len(pvalues)
    order = sorted(range(n), key=lambda i: pvalues[i])
    q = [0.0] * n
    prev = 1.0
    for rank, i in reversed(list(enumerate(order, 1))):
        prev = min(prev, pvalues[i] * n / rank)
        q[i] = prev
    return q


for kind in ['tree', 'BBAA', 'Dmin']:
    path = f'{prefix}_{kind}.txt'
    if not os.path.exists(path):
        continue
    rows = list(csv.DictReader(open(path), delimiter='\t'))
    n = len(rows)
    p = [float(r['p-value']) for r in rows]
    q = bh(p)
    sig = [r for r, qq in zip(rows, q) if qq < alpha]
    n_bonf = sum(x < alpha / n for x in p)

    print(f'\n==== {path} ====')
    print(f'Trios: {n}')
    print(f'p < {alpha} (uncorrected): {sum(x < alpha for x in p)}')
    print(f'Benjamini-Hochberg q < {alpha}: {len(sig)}')
    print(f'Bonferroni p < {alpha / n:.2e}: {n_bonf}')

    print('\nSignificant trios (BH) by populations of P2 and P3:')
    for (a, b), c in Counter(tuple(sorted((pop(r['P2']), pop(r['P3'])))) for r in sig).most_common(15):
        print(f'  {a:8} {b:8} {c}')

    print('\nTop 10 trios:')
    print('  P1\tP2\tP3\tD\tZ\tp\tf4-ratio')
    for r in sorted(sig, key=lambda r: float(r['p-value']))[:10]:
        print('  ' + '\t'.join([r['P1'], r['P2'], r['P3'], r['Dstatistic'], r['Z-score'], r['p-value'], r['f4-ratio']]))

    out = f'{prefix}_{kind}_BH.txt'
    with open(out, 'w') as f:
        f.write('\t'.join(list(rows[0].keys()) + ['q_BH']) + '\n')
        for r, qq in sorted(zip(rows, q), key=lambda x: x[1]):
            f.write('\t'.join(list(r.values()) + [f'{qq:.4g}']) + '\n')
    print(f'\nWritten: {out}')
