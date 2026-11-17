#!/usr/bin/env python3
"""
Make the Dsuite SETS file and guide tree from a popmap and an SVDQuartets tree
(individuals as tips), as produced by add_support_and_pops.py.

SETS file (<VCF sample ID> <tab> <set name>):
  --level individual : each individual is its own set, named <Population>_<ID>
  --level population : set = population
  The outgroup population is named "Outgroup".

Tree (only with --level individual):
  - support values and branch lengths removed;
  - tips renamed <id-prefix><ID> -> <Population>_<ID> (same names as the SETS file);
  - the outgroup individuals replaced by one "Outgroup" tip: (ingroup,Outgroup);
  The outgroup individuals must form one side of the root (tree rooted on the
  outgroup, as saved by run_SVDQuartets.sh).

Usage:
    python3 make_dsuite_inputs.py --popmap popmap.txt --outgroup-pop OUT \\
        --tree results/SVDQ_my_data_besttree_support.tre --out my_SETS
    → my_SETS.txt and my_SETS_tree.nwk

    python3 make_dsuite_inputs.py --popmap popmap.txt --outgroup-pop OUT \\
        --level population --out my_SETS_pops
    → my_SETS_pops.txt (a population tree must be made separately)
"""
import argparse
import re
import sys

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('--popmap', required=True, help='<Population> <ID> file')
parser.add_argument('--outgroup-pop', required=True, help='population name of the outgroup')
parser.add_argument('--tree', help='Newick tree with individuals as tips (needed for --level individual)')
parser.add_argument('--level', choices=['individual', 'population'], default='individual')
parser.add_argument('--id-prefix', default='Ind_', help='prefix of the tip names in the tree (default: Ind_)')
parser.add_argument('--out', required=True, help='output prefix')
args = parser.parse_args()

# ---- popmap (split() also removes Windows line endings) ----
popmap = []
for line in open(args.popmap):
    if line.strip():
        pop, ind = line.split()
        popmap.append((pop, ind))
pops = dict((ind, pop) for pop, ind in popmap)
if args.outgroup_pop not in pops.values():
    sys.exit(f'Outgroup population "{args.outgroup_pop}" not in the popmap')


def set_name(ind):
    pop = pops[ind]
    if pop == args.outgroup_pop:
        return 'Outgroup'
    return f'{pop}_{ind}' if args.level == 'individual' else pop


with open(f'{args.out}.txt', 'w') as f:
    for pop, ind in popmap:
        f.write(f'{ind}\t{set_name(ind)}\n')
print(f'Written: {args.out}.txt ({len(popmap)} samples, '
      f'{len({set_name(i) for _, i in popmap}) - 1} sets + Outgroup)')

if args.level == 'population':
    sys.exit(0)
if not args.tree:
    sys.exit('--tree is needed with --level individual')

# ---- tree ----
newick = open(args.tree).read().strip()
newick = re.sub(r'\)[^(),;]+', ')', newick)        # remove support values / node labels
newick = re.sub(r':[0-9.eE+-]+', '', newick)        # remove branch lengths


def parse(s, i=0):
    if s[i] == '(':
        children, i = [], i + 1
        while True:
            child, i = parse(s, i)
            children.append(child)
            if s[i] == ',':
                i += 1
                continue
            return children, i + 1
    m = re.match(r'[^,();]+', s[i:])
    return m.group(0).strip(), i + m.end()


def tips(node):
    return [node] if isinstance(node, str) else [t for c in node for t in tips(c)]


def to_newick(node):
    if isinstance(node, str):
        ind = node[len(args.id_prefix):] if node.startswith(args.id_prefix) else node
        if ind not in pops:
            sys.exit(f'Tip {node} not in the popmap')
        return set_name(ind)
    return '(' + ','.join(to_newick(c) for c in node) + ')'


tree, _ = parse(newick.rstrip(';'))
outgroup = {args.id_prefix + i for p, i in popmap if p == args.outgroup_pop}
ingroup = [c for c in tree if not set(tips(c)) & outgroup]
outside = [c for c in tree if set(tips(c)) & outgroup]
if len(ingroup) != 1 or any(set(tips(c)) - outgroup for c in outside):
    sys.exit('The outgroup is not on one side of the root: root the tree on the outgroup first')

dsuite_tree = '(' + to_newick(ingroup[0]) + ',Outgroup);'
with open(f'{args.out}_tree.nwk', 'w') as f:
    f.write(dsuite_tree + '\n')

tree_names = set(re.findall(r'[^(),;]+', dsuite_tree))
sets_names = {set_name(i) for _, i in popmap}
if tree_names != sets_names:
    sys.exit(f'Tree and SETS names differ: {sorted(tree_names ^ sets_names)[:10]}')
print(f'Written: {args.out}_tree.nwk ({len(tree_names)} tips, same names as the SETS file)')
