#!/usr/bin/env python3
"""
FEEMS: fast estimation of effective migration surfaces, from a VCF.

Steps (choose with --step):
  graph  build the grid inside the habitat polygon and plot it with the samples.
         Check this plot first: the grid must be fine enough that different
         sampling sites fall on different nodes.
  cv     cross-validation over a range of lambda (smoothing) values; plots the
         CV error and writes the best lambda.
  fit    fit FEEMS with --lamb (or the best lambda from the cv step) and plot
         the migration surface.
  all    graph + cv + fit (default)

Input:
  --vcf     VCF of biallelic SNPs (missing genotypes are replaced by the mean)
  --coords  <sample ID> <tab> <longitude> <tab> <latitude>, no header
  --outer   habitat polygon: longitude latitude, one vertex per line
  --grid    triangular grid shapefile (FEEMS grid_100.shp, or make_feems_grid.R)

Output (in --outdir, prefix --prefix):
  <prefix>_graph.pdf          grid, observed nodes and samples
  <prefix>_cv.pdf / _cv.txt   cross-validation error per lambda
  <prefix>_fit_lamb<L>.pdf    effective migration surface

Usage:
  python3 run_feems.py --vcf my_data.vcf --coords coords_individuals.txt \\
      --outer feems/my_data.outer --grid feems/grid_res8.shp --outdir feems --prefix my_data
"""
import argparse
import gzip
import os
import sys

import numpy as np

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('--vcf', required=True)
parser.add_argument('--coords', required=True, help='<ID> <lon> <lat>, tab-separated, no header')
parser.add_argument('--outer', required=True, help='habitat polygon, lon lat')
parser.add_argument('--grid', required=True, help='grid shapefile (.shp)')
parser.add_argument('--outdir', default='feems')
parser.add_argument('--prefix', default='feems')
parser.add_argument('--step', choices=['graph', 'cv', 'fit', 'all'], default='all')
parser.add_argument('--lamb', type=float, help='lambda for the fit step (default: best lambda from the cv step)')
parser.add_argument('--lamb-min', type=float, default=1e-3, help='smallest lambda tested in cv (default 1e-3)')
parser.add_argument('--lamb-max', type=float, default=1e2, help='largest lambda tested in cv (default 100)')
parser.add_argument('--n-lamb', type=int, default=10, help='number of lambda values in cv (default 10)')
parser.add_argument('--folds', type=int, default=5, help='cross-validation folds (default 5; 0 = leave-one-node-out)')
parser.add_argument('--no-map', action='store_true', help='do not draw coastlines (cartopy downloads them the first time)')
args = parser.parse_args()

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from feems import SpatialGraph, Viz
from feems.utils import prepare_graph_inputs
from feems.cross_validation import run_cv

os.makedirs(args.outdir, exist_ok=True)
out = lambda x: os.path.join(args.outdir, args.prefix + x)

# ---- genotypes from the VCF: copies of the alternative allele (0, 1, 2) ----
opener = gzip.open if args.vcf.endswith('.gz') else open
samples, rows = None, []
with opener(args.vcf, 'rt') as f:
    for line in f:
        if line.startswith('##'):
            continue
        fields = line.rstrip('\n\r').split('\t')
        if line.startswith('#CHROM'):
            samples = fields[9:]
            continue
        if ',' in fields[4]:
            continue                                  # not biallelic
        gt_i = fields[8].split(':').index('GT')
        row = []
        for s in fields[9:]:
            gt = s.split(':')[gt_i].replace('|', '/').split('/')
            row.append(np.nan if (len(gt) != 2 or '.' in gt) else int(gt[0] != '0') + int(gt[1] != '0'))
        rows.append(row)
G = np.array(rows, dtype=float).T                     # individuals x SNPs
# replace missing genotypes by the SNP mean (as in the FEEMS tutorials)
col_mean = np.nanmean(G, axis=0)
idx = np.where(np.isnan(G))
G[idx] = np.take(col_mean, idx[1])
keep = np.nanvar(G, axis=0) > 0                       # drop SNPs that are invariant in these samples
genotypes = G[:, keep]
print(f'n_samples={genotypes.shape[0]}, n_snps={genotypes.shape[1]}')

# ---- coordinates in the VCF order ----
co = {}
with open(args.coords) as f:
    for line in f:
        p = line.split()
        if len(p) >= 3:
            co[p[0]] = (float(p[1]), float(p[2]))
missing = [s for s in samples if s not in co]
if missing:
    sys.exit('Samples without coordinates: ' + ', '.join(missing))
coord = np.array([co[s] for s in samples])
if np.any(np.abs(coord[:, 1]) > 90):
    sys.exit('Latitude > 90: the coordinates file must be ID, longitude, latitude')
outer = np.loadtxt(args.outer)[:, :2]

# ---- graph ----
outer, edges, grid, _ = prepare_graph_inputs(coord=coord, ggrid=args.grid, translated=False,
                                             buffer=0, outer=outer)
sp_graph = SpatialGraph(genotypes, coord, grid, edges, scale_snps=True)
n_sites = len({tuple(c) for c in coord})
print(f'Grid nodes: {grid.shape[0]}, observed nodes: {sp_graph.n_observed_nodes}, sampling sites: {n_sites}')
if sp_graph.n_observed_nodes < n_sites:
    print('Note: several sampling sites share a node. If they are different populations, '
          'use a finer grid (make_feems_grid.R with a higher RES).')

center = coord.mean(axis=0)
projection = ccrs.EquidistantConic(central_longitude=center[0], central_latitude=center[1])


def new_viz(title):
    fig = plt.figure(dpi=300)
    ax = fig.add_subplot(1, 1, 1, projection=projection)
    ax.set_title(title, fontsize=8)
    v = Viz(ax, sp_graph, projection=projection, edge_width=0.5, edge_alpha=1, edge_zorder=100,
            sample_pt_size=10, obs_node_size=7.5, sample_pt_color='black', cbar_font_size=8)
    if not args.no_map:
        v.draw_map()
    return fig, v


if args.step in ('graph', 'all'):
    fig, v = new_viz(f'{grid.shape[0]} nodes, {sp_graph.n_observed_nodes} observed')
    v.draw_samples()
    v.draw_edges(use_weights=False)
    v.draw_obs_nodes(use_ids=False)
    fig.savefig(out('_graph.pdf'))
    print('Graph plot:', out('_graph.pdf'))

# ---- cross-validation ----
cv_file = out('_cv.txt')
if args.step in ('cv', 'all'):
    lamb_grid = np.geomspace(args.lamb_min, args.lamb_max, args.n_lamb, endpoint=True)[::-1]   # large to small: warm start
    n_folds = None if args.folds == 0 else min(args.folds, sp_graph.n_observed_nodes)
    cv_err = run_cv(sp_graph, lamb_grid, n_folds=n_folds, factr=1e10, outer_verbose=False)
    mean_cv = np.mean(cv_err, axis=0)[:, 0]
    best = float(lamb_grid[np.argmin(mean_cv)])
    np.savetxt(cv_file, np.column_stack([lamb_grid, mean_cv]), header='lambda\tmean_cv_error', fmt='%.6g',
               delimiter='\t', comments='')
    fig, ax = plt.subplots(dpi=200)
    ax.plot(lamb_grid, mean_cv, '-o')
    ax.set_xscale('log'); ax.set_xlabel('lambda'); ax.set_ylabel('CV error')
    ax.axvline(best, color='orange')
    fig.savefig(out('_cv.pdf'))
    print(f'Best lambda (lowest CV error): {best:.4g}  ->', cv_file)
    if best in (lamb_grid.min(), lamb_grid.max()):
        print('Note: the best lambda is at the edge of the tested range; widen --lamb-min/--lamb-max.')

# ---- fit ----
if args.step in ('fit', 'all'):
    lamb = args.lamb
    if lamb is None:
        if not os.path.exists(cv_file):
            sys.exit('Give --lamb or run the cv step first')
        tab = np.loadtxt(cv_file, skiprows=1)
        lamb = float(tab[np.argmin(tab[:, 1]), 0])
    sp_graph.fit(lamb=lamb, optimize_q=None)   # one residual variance for all nodes (as in run_cv)
    fig, v = new_viz(f'lambda = {lamb:.4g}')
    v.draw_edges(use_weights=True)
    v.draw_obs_nodes(use_ids=False)
    v.draw_samples()
    v.draw_edge_colorbar()
    fig.savefig(out(f'_fit_lamb{lamb:.4g}.pdf'))
    print('Migration surface:', out(f'_fit_lamb{lamb:.4g}.pdf'))
