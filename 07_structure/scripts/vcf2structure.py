#!/usr/bin/env python3
"""
Convert a VCF of biallelic SNPs to the input files of STRUCTURE and
fastSTRUCTURE (run with structure_threader), with the individuals ordered by
population.

Output (<out> = --out prefix):
  <out>.str            STRUCTURE format: a row of marker names, then two rows
                       per individual: <ID> <pop number> <genotypes...>
                       (LABEL=1, POPDATA=1, MARKERNAMES=1, EXTRACOLS=0, ONEROWPERIND=0)
  <out>_fast.str       fastSTRUCTURE format: two rows per individual, six
                       leading columns (ignored by fastSTRUCTURE), no marker names
  <out>_popfile.txt    <population> <number of individuals> <order>, for
                       structure_threader --pop (labels of the plots)
  <out>_individuals.txt <ID> <tab> <population>, in the order of the .str files
                       (labels for plot_structure.R)

Alleles are coded 0 (reference) and 1 (alternative); missing = -9.

Usage:
    python3 vcf2structure.py --vcf my_data.vcf --popmap popmap.txt --out my_data
"""
import argparse
import gzip
import sys

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('--vcf', required=True, help='VCF (plain or .gz) with biallelic SNPs')
parser.add_argument('--popmap', required=True, help='<population> <tab> <sample ID>')
parser.add_argument('--out', required=True, help='output prefix')
args = parser.parse_args()

# ---- popmap: population order = order of first appearance ----
pop_of = {}
pop_order = []
with open(args.popmap) as f:
    for line in f:
        parts = line.split()          # split() also removes Windows line endings
        if len(parts) < 2:
            continue
        pop, ind = parts[0], parts[1]
        if ind in pop_of:
            sys.exit(f'Sample {ind} appears twice in the popmap')
        pop_of[ind] = pop
        if pop not in pop_order:
            pop_order.append(pop)

# ---- read VCF ----
opener = gzip.open if args.vcf.endswith('.gz') else open
samples = None
markers = []
genos = []                            # one list per SNP: (allele1, allele2) per sample
skipped = 0
with opener(args.vcf, 'rt') as f:
    for line in f:
        if line.startswith('##'):
            continue
        fields = line.rstrip('\n\r').split('\t')
        if line.startswith('#CHROM'):
            samples = fields[9:]
            continue
        if samples is None:
            sys.exit('No #CHROM header line in the VCF')
        alt = fields[4]
        if ',' in alt or len(fields[3]) != 1 or len(alt) != 1:
            skipped += 1              # not a biallelic SNP
            continue
        fmt = fields[8].split(':')
        if 'GT' not in fmt:
            sys.exit(f'No GT field at {fields[0]}:{fields[1]}')
        gt_i = fmt.index('GT')
        row = []
        for s in fields[9:]:
            gt = s.split(':')[gt_i].replace('|', '/').split('/')
            if len(gt) != 2 or '.' in gt:
                row.append(('-9', '-9'))
            else:
                row.append((gt[0], gt[1]))
        genos.append(row)
        markers.append(f'{fields[0]}_{fields[1]}')

if samples is None:
    sys.exit('No #CHROM header line in the VCF')
missing_pop = [s for s in samples if s not in pop_of]
if missing_pop:
    sys.exit('Samples missing from the popmap: ' + ', '.join(missing_pop))
not_in_vcf = [s for s in pop_of if s not in samples]
if not_in_vcf:
    print(f'Note: {len(not_in_vcf)} samples of the popmap are not in the VCF (ignored)', file=sys.stderr)

# ---- order individuals by population ----
pop_number = {p: i + 1 for i, p in enumerate(pop_order)}
idx = sorted(range(len(samples)), key=lambda i: (pop_number[pop_of[samples[i]]], i))

with open(f'{args.out}.str', 'w') as st, open(f'{args.out}_fast.str', 'w') as fs:
    st.write('\t'.join(markers) + '\n')
    for i in idx:
        ind = samples[i]
        p = pop_number[pop_of[ind]]
        for allele in (0, 1):
            g = [genos[m][i][allele] for m in range(len(markers))]
            st.write(f'{ind}\t{p}\t' + '\t'.join(g) + '\n')
            fs.write(f'{ind}\t{p}\t0\t0\t0\t0\t' + '\t'.join(g) + '\n')

used_pops = [p for p in pop_order if any(pop_of[s] == p for s in samples)]
with open(f'{args.out}_popfile.txt', 'w') as f:
    for p in used_pops:
        n = sum(1 for s in samples if pop_of[s] == p)
        f.write(f'{p}\t{n}\t{pop_number[p]}\n')

with open(f'{args.out}_individuals.txt', 'w') as f:
    for i in idx:
        f.write(f'{samples[i]}\t{pop_of[samples[i]]}\n')

print(f'Individuals (NUMINDS): {len(samples)}')
print(f'Loci (NUMLOCI):        {len(markers)}')
if skipped:
    print(f'Skipped {skipped} sites that are not biallelic SNPs')
print(f'Populations:           {len(used_pops)}')
print(f'Written: {args.out}.str, {args.out}_fast.str, {args.out}_popfile.txt, {args.out}_individuals.txt')
