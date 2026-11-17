# FEEMS: effective migration surface with cross-validation.
#
# Run from the analysis folder (10_feems/):
#   conda activate feems
#   python scripts/run_feems.py
# (or open it in Jupyter / VS Code and run it block by block)

import numpy as np
import matplotlib
matplotlib.use("Agg")                  # save the plots to files
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from feems import SpatialGraph, Viz
from feems.utils import prepare_graph_inputs
from feems.cross_validation import run_cv

# ---- settings: edit these lines ----
vcf_file = "data/my_data-miss0.9-maf0.01-th100.recode.vcf"
coords_file = "data/coords.txt"         # <sample ID> <longitude> <latitude>
outer_file = "data/my_data.outer"       # polygon of the study area (longitude latitude)
grid_file = "data/grid_res8.shp"        # grid made with make_feems_grid.R
draw_map = True                         # False if cartopy cannot download the coastlines (no internet)

# ---- 1. read the genotypes from the VCF (0, 1, 2 copies of the alternative allele) ----
samples = []
rows = []
for line in open(vcf_file):
    if line.startswith("##"):
        continue
    fields = line.strip().split("\t")
    if line.startswith("#CHROM"):
        samples = fields[9:]
        continue
    row = []
    for sample_field in fields[9:]:
        gt = sample_field.split(":")[0].replace("|", "/")
        if gt == "0/0":
            row.append(0)
        elif gt in ("0/1", "1/0"):
            row.append(1)
        elif gt == "1/1":
            row.append(2)
        else:
            row.append(np.nan)          # missing
    rows.append(row)

genotypes = np.array(rows, dtype=float).T          # individuals x SNPs
# replace missing genotypes by the mean of the SNP (as in the FEEMS tutorials)
snp_means = np.nanmean(genotypes, axis=0)
missing = np.where(np.isnan(genotypes))
genotypes[missing] = snp_means[missing[1]]
print("individuals:", genotypes.shape[0], " SNPs:", genotypes.shape[1])

# ---- 2. coordinates in the same order as the VCF ----
coords_table = {}
for line in open(coords_file):
    sample_id, lon, lat = line.split()[:3]
    coords_table[sample_id] = (float(lon), float(lat))
coord = np.array([coords_table[s] for s in samples])   # error here = a sample without coordinates
outer = np.loadtxt(outer_file)

# ---- 3. build the graph and plot it ----
outer, edges, grid, _ = prepare_graph_inputs(coord=coord, ggrid=grid_file, translated=False,
                                             buffer=0, outer=outer)
sp_graph = SpatialGraph(genotypes, coord, grid, edges, scale_snps=True)
print("nodes in the grid:", grid.shape[0], " nodes with samples:", sp_graph.n_observed_nodes)

projection = ccrs.EquidistantConic(central_longitude=coord[:, 0].mean(),
                                   central_latitude=coord[:, 1].mean())

fig = plt.figure(dpi=300)
ax = fig.add_subplot(1, 1, 1, projection=projection)
v = Viz(ax, sp_graph, projection=projection, edge_width=0.5, edge_alpha=1, edge_zorder=100,
        sample_pt_size=10, obs_node_size=7.5, sample_pt_color="black", cbar_font_size=8)
if draw_map:
    v.draw_map()
v.draw_samples()
v.draw_edges(use_weights=False)
v.draw_obs_nodes(use_ids=False)
fig.savefig("results/feems_graph.pdf")

# ---- 4. cross-validation: choose lambda (the smoothing) ----
lamb_grid = np.geomspace(1e-3, 1e2, 10)[::-1]    # 10 values from 100 to 0.001
cv_err = run_cv(sp_graph, lamb_grid, n_folds=5, factr=1e10)
mean_cv_err = np.mean(cv_err, axis=0)[:, 0]
best_lamb = float(lamb_grid[np.argmin(mean_cv_err)])
print("lambda and CV error:")
for lamb, err in zip(lamb_grid, mean_cv_err):
    print(f"  {lamb:10.4f}  {err:.6f}")
print("best lambda:", best_lamb)

fig, ax = plt.subplots(dpi=200)
ax.plot(lamb_grid, mean_cv_err, "-o")
ax.set_xscale("log")
ax.set_xlabel("lambda")
ax.set_ylabel("CV error")
ax.axvline(best_lamb, color="orange")
fig.savefig("results/feems_cv.pdf")

# ---- 5. fit FEEMS with the best lambda and plot the migration surface ----
sp_graph.fit(lamb=best_lamb, optimize_q=None)

fig = plt.figure(dpi=300)
ax = fig.add_subplot(1, 1, 1, projection=projection)
v = Viz(ax, sp_graph, projection=projection, edge_width=0.5, edge_alpha=1, edge_zorder=100,
        sample_pt_size=10, obs_node_size=7.5, sample_pt_color="black", cbar_font_size=8)
if draw_map:
    v.draw_map()
v.draw_edges(use_weights=True)
v.draw_obs_nodes(use_ids=False)
v.draw_samples()
v.draw_edge_colorbar()
fig.savefig("results/feems_migration.pdf")
print("Done: results/feems_graph.pdf, results/feems_cv.pdf, results/feems_migration.pdf")
