#!/usr/bin/env python3
"""
Check whether two groups of populations (or every pair of populations) have
fixed differences, before running ancestry painting.

Uses biallelic SNPs only (single-base REF and ALT). Reads plain or gzipped VCF.

1) Two groups (e.g. the two putative parents), defined by regular expressions
   on the population names of the popmap:
       python3 count_fixed_differences.py my_data.vcf popmap.txt --p1 '^PopA' --p2 '^PopB'
   → distribution of the allele-frequency differences between the two groups.

2) All pairs of populations:
       python3 count_fixed_differences.py my_data.vcf popmap.txt --pairs
   → number of fixed differences for every pair of populations.

A site is used for a group/population when at least --min-data (default 0.5)
of its alleles are genotyped. "Fixed" = one side has only the reference
allele and the other only the alternative allele.
"""
import argparse
import gzip
import itertools
import re
import sys

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('vcf')
parser.add_argument('popmap', help='<Population> <ID> file')
parser.add_argument('--p1', help='regular expression for the populations of group 1')
parser.add_argument('--p2', help='regular expression for the populations of group 2')
parser.add_argument('--pairs', action='store_true', help='count fixed differences for all pairs of populations')
parser.add_argument('--min-data', type=float, default=0.5, help='min. proportion of alleles genotyped (default 0.5)')
args = parser.parse_args()
if not args.pairs and not (args.p1 and args.p2):
    sys.exit('Use --p1 and --p2, or --pairs')

pops = {}
for line in open(args.popmap):
    if line.strip():
        pop, ind = line.split()
        pops[ind] = pop

if args.pairs:
    groups = {p: p for p in sorted(set(pops.values()))}
    groups = {p: (lambda name, p=p: name == p) for p in groups}
else:
    groups = {'group1': lambda name: re.search(args.p1, name),
              'group2': lambda name: re.search(args.p2, name)}

opener = gzip.open if args.vcf.endswith('.gz') else open
idx = {}
freqs_per_site = []
n_snps = 0
with opener(args.vcf, 'rt') as f:
    for line in f:
        if line.startswith('##'):
            continue
        fields = line.rstrip('\n').split('\t')
        if line.startswith('#'):
            for g, test in groups.items():
                idx[g] = [i for i, s in enumerate(fields) if s in pops and test(pops[s])]
                if not idx[g]:
                    sys.exit(f'No samples for {g}')
            for g in idx:
                print(f'{g}: {len(idx[g])} samples', file=sys.stderr)
            continue
        if len(fields[3]) != 1 or len(fields[4]) != 1:
            continue
        n_snps += 1
        gt = fields[8].split(':').index('GT')
        freqs = {}
        for g, cols in idx.items():
            alleles = [a for i in cols for a in fields[i].split(':')[gt].replace('|', '/').split('/') if a not in ('.', '')]
            if len(alleles) >= args.min_data * 2 * len(cols):
                freqs[g] = alleles.count('1') / len(alleles)
        freqs_per_site.append(freqs)

print(f'Biallelic SNPs: {n_snps}')

if not args.pairs:
    diffs = [abs(f['group1'] - f['group2']) for f in freqs_per_site if len(f) == 2]
    print(f'SNPs with data in both groups: {len(diffs)}')
    if diffs:
        print(f'Max. allele-frequency difference: {max(diffs):.3f}')
        for t in [1.0, 0.95, 0.9, 0.8, 0.6]:
            print(f'  difference >= {t}: {sum(d >= t - 1e-9 for d in diffs)}')
else:
    print('pop1\tpop2\tfixed\tsites_compared')
    rows = []
    for a, b in itertools.combinations(sorted(groups), 2):
        both = [f for f in freqs_per_site if a in f and b in f]
        fixed = sum({f[a], f[b]} == {0.0, 1.0} for f in both)
        rows.append((fixed, a, b, len(both)))
    for fixed, a, b, n in sorted(rows, reverse=True):
        print(f'{a}\t{b}\t{fixed}\t{n}')
