#!/usr/bin/env python3
"""
Remove the columns of a PHYLIP alignment that RAxML-ng considers invariant,
so that the ascertainment bias correction (+ASC_LEWIS) can be used.

A column is invariant when one nucleotide is compatible with all genotypes of
the column. With heterozygotes written as IUPAC codes, a SNP column with only
A and R (A/G) is invariant: every individual can carry an A.

Usage:
    python3 remove_invariant_sites.py input.phy output.phy
"""
import sys

if len(sys.argv) != 3:
    sys.exit(__doc__)

# nucleotides that each IUPAC code can be
IUPAC = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
         'R': 'AG', 'Y': 'CT', 'S': 'CG', 'W': 'AT', 'K': 'GT', 'M': 'AC',
         'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG'}

# read the PHYLIP file (sequential format: one line per individual)
with open(sys.argv[1]) as f:
    n_ind, n_sites = map(int, f.readline().split())
    names, seqs = [], []
    for line in f:
        if line.strip():
            name, seq = line.split()
            names.append(name)
            seqs.append(seq.upper())

keep = []
for i in range(n_sites):
    common = set('ACGT')
    for seq in seqs:
        base = seq[i]
        if base in IUPAC:                     # N, - and ? (missing) are ignored
            common &= set(IUPAC[base])
    if not common:                            # no nucleotide shared by all: variable column
        keep.append(i)

with open(sys.argv[2], 'w') as out:
    out.write(f'{n_ind} {len(keep)}\n')
    width = max(len(n) for n in names) + 3
    for name, seq in zip(names, seqs):
        out.write(name.ljust(width) + ''.join(seq[i] for i in keep) + '\n')

print(f'{n_sites} sites in, {n_sites - len(keep)} invariant removed, {len(keep)} kept: {sys.argv[2]}')
