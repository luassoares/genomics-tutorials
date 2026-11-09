#!/usr/bin/env python3
"""
Add bootstrap support values and population names to the SVDQuartets trees.

PAUP* does not write the support values with savetrees for the SVDQuartets
consensus tree, so they are calculated here from the bootstrap trees
(<prefix>_bootstrap.tre): support of a node = % of bootstrap trees that
contain the same bipartition.

Usage:
    python3 add_support_and_pops.py <results_dir/prefix> <popmap> [--id-prefix Ind_]

Example:
    python3 add_support_and_pops.py results/SVDQ_my_data ../popmap.txt

Input (from run_SVDQuartets.sh):
    <prefix>_bootstrap.tre   bootstrap trees
    <prefix>_consensus.tre   majority-rule consensus
    <prefix>_besttree.tre    best (fully resolved) tree

Output (Newick, support as node labels):
    <prefix>_consensus_support.tre        tips = <id-prefix><ID>
    <prefix>_consensus_support_pops.tre   tips = <Population>_<ID>
    <prefix>_besttree_support.tre
    <prefix>_besttree_support_pops.tre

Popmap: two columns, tab/space separated: <Population> <ID>
(IDs as in the VCF, i.e. without the prefix added to the NEXUS taxon names;
use --id-prefix "" if no prefix was added).
"""
import argparse
import os
import re
import sys
from collections import Counter

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('prefix', help='results prefix, e.g. results/SVDQ_my_data')
parser.add_argument('popmap', help='<Population> <ID> file')
parser.add_argument('--id-prefix', default='Ind_', help='prefix of the taxon names in the trees (default: Ind_)')
args = parser.parse_args()
prefix, popmap, id_prefix = args.prefix, args.popmap, args.id_prefix

# ---- popmap: <id_prefix><ID> -> <Population>_<ID> ----
pops = {}
for line in open(popmap):
    if line.strip():
        pop, ind = line.split()
        pops[id_prefix + ind] = f'{pop}_{ind}'


def read_nexus_trees(path):
    """Return the trees of a PAUP* NEXUS tree file, with the tip numbers replaced by taxon names."""
    txt = open(path).read()
    names = {}
    taxlabels = re.search(r'Taxlabels(.*?);', txt, re.S | re.I)
    if taxlabels:
        names = {str(k + 1): v for k, v in enumerate(taxlabels.group(1).split())}
    translate = re.search(r'Translate(.*?);', txt, re.S | re.I)
    if translate and translate.group(1).strip():
        names.update({k: v.strip("'") for k, v in (e.split() for e in translate.group(1).split(','))})
    trees = re.findall(r'tree\s+\S+\s*=\s*(?:\[[^\]]*\]\s*)*(\(.*?\));', txt, re.I)
    return [re.sub(r'(?<=[(,])\s*([^(),;:\s]+)', lambda m: names.get(m.group(1), m.group(1)), t)
            for t in trees]


def tokens(newick):
    return re.findall(r'\(|\)[^(),;]*|[^(),;]+|,|;', newick)


def clades(newick):
    """All clades (lists of tip names) of a Newick tree."""
    out, stack = [], [[]]
    for t in tokens(newick):
        if t == '(':
            stack.append([])
        elif t[0] == ')':
            clade = stack.pop()
            out.append(clade)
            stack[-1].extend(clade)
        elif t not in (',', ';'):
            stack[-1].append(t.strip())
    return out


# ---- count bipartitions in the bootstrap trees ----
boots = read_nexus_trees(f'{prefix}_bootstrap.tre')
all_tips = {tip for c in clades(boots[0]) for tip in c}
ref = sorted(all_tips)[0]


def bipartition(clade):
    """Unrooted bipartition, always written as the side without the reference tip."""
    s = frozenset(clade)
    return s if ref not in s else frozenset(all_tips - s)


counts = Counter()
for tree in boots:
    counts.update({bipartition(c) for c in clades(tree)})
n = len(boots)
print(f'{n} bootstrap trees read from {prefix}_bootstrap.tre')


def add_support(newick):
    out, stack = [], [[]]
    for t in tokens(newick):
        if t == '(':
            stack.append([])
            out.append('(')
        elif t[0] == ')':
            clade = stack.pop()
            stack[-1].extend(clade)
            is_root = len(stack) == 1
            out.append(')' + ('' if is_root else str(round(100 * counts[bipartition(clade)] / n))))
        elif t in (',', ';'):
            out.append(t)
        else:
            stack[-1].append(t.strip())
            out.append(t.strip())
    return ''.join(out)


# ---- write the trees ----
for kind in ['consensus', 'besttree']:
    path = f'{prefix}_{kind}.tre'
    if not os.path.exists(path):
        print(f'{path} not found, skipped')
        continue
    tree = add_support(read_nexus_trees(path)[0])
    with open(f'{prefix}_{kind}_support.tre', 'w') as f:
        f.write(tree + '\n')
    tips = {t for c in clades(tree) for t in c}
    missing = tips - set(pops)
    if missing:
        sys.exit(f'Not in popmap: {" ".join(sorted(missing))}')
    renamed = re.sub(r'(?<=[(,])([^(),;]+)', lambda m: pops[m.group(1)], tree)
    with open(f'{prefix}_{kind}_support_pops.tre', 'w') as f:
        f.write(renamed + '\n')
    print(f'Written: {prefix}_{kind}_support.tre and {prefix}_{kind}_support_pops.tre')
