#!/usr/bin/env python3
"""
Choose the individuals for Dfoil: in each population, the N individuals with
the least missing data.

Usage:
    python3 pick_individuals.py popmap.txt my_data.imiss N > names.txt

popmap.txt     <population> <tab> <sample ID>
my_data.imiss  from: vcftools --vcf my_data.vcf --missing-indv --out my_data
N              individuals per population (e.g. 3)

The outgroup populations (names starting with "Outgroup") are not written:
the outgroup goes in a separate file (see the tutorial).
A table with the chosen individuals and their missing data is printed to
the screen (stderr).
"""
import sys
from collections import defaultdict

if len(sys.argv) != 4:
    sys.exit(__doc__)

popmap_file, imiss_file, n = sys.argv[1], sys.argv[2], int(sys.argv[3])

# population of each sample
pop_of = {}
with open(popmap_file) as f:
    for line in f:
        line = line.strip()
        if line:
            pop, sample = line.split()[:2]
            pop_of[sample] = pop

# missing data of each sample (column F_MISS of vcftools --missing-indv)
fmiss = {}
with open(imiss_file) as f:
    header = f.readline().split()
    col = header.index("F_MISS")
    for line in f:
        fields = line.split()
        fmiss[fields[0]] = float(fields[col])

by_pop = defaultdict(list)
for sample, pop in pop_of.items():
    if sample not in fmiss:
        print(f"Warning: {sample} is in the popmap but not in the VCF", file=sys.stderr)
        continue
    by_pop[pop].append((fmiss[sample], sample))

print("population\tsample\tF_MISS", file=sys.stderr)
total = 0
for pop in sorted(by_pop):
    if pop.startswith("Outgroup"):
        continue
    for miss, sample in sorted(by_pop[pop])[:n]:
        print(sample)
        print(f"{pop}\t{sample}\t{miss:.3f}", file=sys.stderr)
        total += 1
print(f"\n{total} individuals chosen", file=sys.stderr)
