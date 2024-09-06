#!/usr/bin/env python
# Author: M. Vochteloo (Adapted from 'https://github.com/sc-eQTLgen-consortium/WG3-pipeline-QTL/blob/master/scripts/testDataset_Base.R' by M.J. Bonder)

import argparse
import gzip
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(description="")
parser.add_argument("--metadata", required=True, type=str, help="")
parser.add_argument("--ancestry", required=True, type=str, help="")
parser.add_argument("--cell_level", required=True, type=str, help="")
parser.add_argument("--cell_type", required=True, type=str, help="")
parser.add_argument("--qc", required=True, type=str, help="")
parser.add_argument("--donor_cell_threshold", required=False, type=float, default=None, help="")
parser.add_argument("--principal_components", required=True, type=str, help="")
parser.add_argument("--counts_threshold", required=False, type=float, default=None, help="")
parser.add_argument("--full_kinship", required=True, type=str, help="")
parser.add_argument("--kinship", required=True, type=str, help="")
parser.add_argument("--kinship_threshold", required=False, type=float, default=None, help="")
parser.add_argument("--smf", required=True, type=str, help="")
parser.add_argument("--individual_aggregate", required=False, default="Assignment", type=str, help="")
parser.add_argument("--sample_aggregate", required=False, default="Assignment_Run_Lane", type=str, help="")
parser.add_argument("--data_out", required=True, type=str, help="The output directory where results will be saved.")
parser.add_argument("--plot_out", required=True, type=str, help="The output directory where plots will be saved.")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")


def gzopen(file, mode="r"):
    if file.endswith(".gz"):
        return gzip.open(file, mode + 't')
    else:
        return open(file, mode)


MEDIAN_COLOR = "#000000"
MAD_THRESHOLDS = {1.: "#1a1a1a", 2.: "#4d4d4d", 3.: "#7f7f7f", 4.: "#b3b3b3", 5.: "#e5e5e5"}
SELECTED_COLOR = "#B22222"

DONOR_CELL_THRESHOLDS = MAD_THRESHOLDS.copy()
if args.donor_cell_threshold is not None:
    if args.donor_cell_threshold in DONOR_CELL_THRESHOLDS:
        DONOR_CELL_THRESHOLDS[args.donor_cell_threshold] = SELECTED_COLOR
    else:
        DONOR_CELL_THRESHOLDS[args.donor_cell_threshold] = SELECTED_COLOR

COUNTS_THRESHOLDS = MAD_THRESHOLDS.copy()
if args.counts_threshold is not None:
    if args.counts_threshold in COUNTS_THRESHOLDS:
        COUNTS_THRESHOLDS[args.counts_threshold] = SELECTED_COLOR
    else:
        COUNTS_THRESHOLDS[args.counts_threshold] = SELECTED_COLOR

KINSHIP_THRESHOLDS = MAD_THRESHOLDS.copy()
if args.kinship_threshold is not None:
    if args.kinship_threshold in KINSHIP_THRESHOLDS:
        KINSHIP_THRESHOLDS[args.kinship_threshold] = SELECTED_COLOR
    else:
        KINSHIP_THRESHOLDS[args.kinship_threshold] = SELECTED_COLOR

individual_aggregate = args.individual_aggregate
sample_aggregate = args.sample_aggregate
if args.individual_aggregate == args.sample_aggregate:
    individual_aggregate = individual_aggregate + "_ind"
    sample_aggregate = sample_aggregate + "_sample"


def calculate_mad_thresholds(data, constant=1.4826, threshold=1.):
    median = np.median(data)
    lower_threshold = median - (np.median(np.abs(data - median))) * constant * threshold
    upper_threshold = median + (np.median(np.abs(data - median))) * constant * threshold
    return lower_threshold, upper_threshold


def filter_on_mad(df, x, y=None, threshold=1.):
    if y is None:
        y = x
    if threshold is None:
        df["include"] = True
        return df

    x_lower_threshold, x_upper_threshold = calculate_mad_thresholds(data=df[x], threshold=threshold)
    y_lower_threshold, y_upper_threshold = calculate_mad_thresholds(data=df[y], threshold=threshold)

    df["include"] = False
    df.loc[(df[x] > x_lower_threshold) & (df[x] < x_upper_threshold) &
           (df[y] > y_lower_threshold) & (df[y] < y_upper_threshold), "include"] = True

    return df


def plot_histogram(ax, y, thresholds=None, title='', xlabel='', ylabel=''):
    if thresholds is None:
        thresholds = MAD_THRESHOLDS

    ax.hist(y, color='grey', linewidth=0, density=True)
    yl = ax.get_ylim()
    ax.set_ylim((yl[0], yl[1] * 1.1))

    add_mad_lines(ax=ax, values=y, thresholds=thresholds, horizontal=False)

    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)


def add_mad_lines(ax, values, thresholds, horizontal=False, alpha=0.4):
    outer_bound = ax.get_xlim()
    line_func = ax.axvline
    if horizontal:
        outer_bound = ax.get_ylim()
        line_func = ax.axhline

    median = np.median(values)

    line_func(median, ls='--', color=MEDIAN_COLOR, alpha=min(1., alpha * 2), zorder=-1)

    for mult, color in thresholds.items():
        lower_threshold, upper_threshold = calculate_mad_thresholds(data=values, threshold=mult)

        alpha_mult = 1
        if color == SELECTED_COLOR:
            alpha_mult = 2

        if upper_threshold < outer_bound[1]:
            line_func(upper_threshold, ls='--', color=color, alpha=min(1., alpha * alpha_mult), zorder=-1)
        if lower_threshold > outer_bound[0]:
            line_func(lower_threshold, ls='--', color=color, alpha=min(1., alpha * alpha_mult), zorder=-1)


def plot_legend(ax, df, x, y=None, thresholds=None):
    if thresholds is None:
        thresholds = MAD_THRESHOLDS

    ax.set_axis_off()

    handles = []
    handles.append((-1, mpatches.Patch(color=MEDIAN_COLOR, label="Median")))
    for threshold, color in thresholds.items():
        tmp_df = filter_on_mad(df=df, x=x, y=y, threshold=threshold)
        n = tmp_df.loc[tmp_df["include"], :].shape[0]
        handles.append((threshold, mpatches.Patch(color=color, label="±{} MAD: N = {:,}".format(threshold, n))))
        del tmp_df

    handles.sort(key=lambda x: x[0])
    handles = [handle[1] for handle in handles]
    ax.legend(handles=handles, loc='center', fontsize=10)


def plot_embedding(ax, embedding, x, y, alpha=1., thresholds=None, z=None, title='', include_colorbar=True):
    if thresholds is None:
        thresholds = MAD_THRESHOLDS

    x_data = embedding.loc[:, x].to_numpy()
    y_data = embedding.loc[:, y].to_numpy()
    xl = (x_data.min() - x_data.ptp() * .05, x_data.max() + x_data.ptp() * 0.05)
    yl = (y_data.min() - y_data.ptp() * .05, y_data.max() + y_data.ptp() * 0.05)

    c = "black"
    vmin = None
    vmax = None
    if z:
        color_dat = embedding.loc[:, z].to_numpy()
        vmin = color_dat.min()
        vmax = color_dat.max()

        o = np.argsort(color_dat)
        x_data = x_data[o]
        y_data = y_data[o]
        c = color_dat[o]

    pp = plot_scatterplot(ax=ax, x_data=x_data, y_data=y_data, c=c, vmin=vmin, vmax=vmax, alpha=alpha)

    add_mad_lines(ax=ax, values=x_data, thresholds=thresholds, horizontal=False)
    add_mad_lines(ax=ax, values=y_data, thresholds=thresholds, horizontal=True)

    ax.set_xlim(xl)
    ax.set_ylim(yl)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(title)
    ax.set_xlabel(x)
    ax.set_ylabel(y)

    if include_colorbar:
        fig.colorbar(pp, ax=ax)


def plot_scatterplot(ax, x_data, y_data, s=30, c="black", vmin=None, vmax=None, alpha=1.):
    return ax.scatter(x_data, y_data,
                      s=s,
                      edgecolors=None,
                      c=c,
                      cmap=darken_cmap(plt.cm.Reds, 0.9),
                      vmin=vmin,
                      vmax=vmax,
                      alpha=alpha)


def darken_cmap(cmap, scale_factor):
    cdat = np.zeros((cmap.N, 4))
    for ii in range(cdat.shape[0]):
        curcol = cmap(ii)
        cdat[ii, 0] = curcol[0] * scale_factor
        cdat[ii, 1] = curcol[1] * scale_factor
        cdat[ii, 2] = curcol[2] * scale_factor
        cdat[ii, 3] = 1
    cmap = cmap.from_list(cmap.N, cdat)
    return cmap


################################################################################

print("Loading SMF file")
smf = pd.read_csv(args.smf, sep="\t", header=None, index_col=None, low_memory=False)
smf.columns = [individual_aggregate, sample_aggregate]

################################################################################

print("Plotting")
plt.rc('font', size=14)
plt.rcParams['pdf.fonttype'] = 42
fig, axs = plt.subplots(3, 3, figsize=(12, 12), dpi=150, gridspec_kw={"width_ratios": [0.45, 0.45, 0.05]})

################################################################################

print("Loading metadata and applying QC filters")
metadata = pd.read_csv(args.metadata, sep="\t", header=0, index_col=None)
# THIS HAS TO BE IDENTICAL TO create_expression_matrices.R filter except for the cell type filter!!!!!!!!!!!
# ps. the fact that individual-sample aggregate is missing here is fine since we only use that to completely remove
# samples and those would not be in this metadata file anymore.
metadata = metadata.loc[(metadata["DropletType"] == "singlet") &
                        (metadata["tag"] == "NotOutlier") &
                        (metadata["Provided_Ancestry"] == args.ancestry) &
                        (metadata["cell_treatment"] == "UT"), :]

# Counts per sample.
count_per_ind = metadata.loc[:, [args.sample_aggregate, "Barcode"]].groupby(args.sample_aggregate).count()
count_per_ind.columns = ["Total"]

# Cell type counts per sample.
ct_count_per_ind = metadata.loc[metadata[args.cell_level] == args.cell_type, [args.sample_aggregate, "Barcode"]].groupby([args.sample_aggregate]).count()
ct_count_per_ind.columns = ["CellCount"]

# Combine.
count_per_ind_df = count_per_ind.merge(ct_count_per_ind, left_index=True, right_index=True)
count_per_ind_df["Fraction"] = count_per_ind_df["CellCount"] / count_per_ind_df["Total"]
del count_per_ind, ct_count_per_ind

plot_histogram(
    ax=axs[0, 0],
    y=count_per_ind_df["Fraction"].to_numpy(),
    thresholds=DONOR_CELL_THRESHOLDS,
    title='Cell fraction'
)
plot_histogram(
    ax=axs[0, 1],
    y=count_per_ind_df["CellCount"].to_numpy(),
    thresholds=DONOR_CELL_THRESHOLDS,
    title='Cell count'
)
plot_legend(
    ax=axs[0, 2],
    df=count_per_ind_df,
    x='CellCount',
    thresholds=DONOR_CELL_THRESHOLDS
)

include_count_per_ind_df = filter_on_mad(df=count_per_ind_df, x="Fraction", threshold=args.donor_cell_threshold)
include_count_per_ind_df = include_count_per_ind_df.loc[:, ["include"]].copy()
include_count_per_ind_df.columns = ["cells"]

smf = smf.merge(include_count_per_ind_df, left_on=sample_aggregate, right_index=True, how="left")

################################################################################

print("Processing Principal Components")
principal_components = pd.read_csv(args.principal_components, sep="\t", header=0, index_col=0)
principal_components = principal_components.merge(count_per_ind_df[["CellCount"]], how="left", left_index=True, right_index=True)
principal_components = filter_on_mad(df=principal_components, x="PC1", y="PC2", threshold=args.counts_threshold)
principal_components["hue"] = principal_components["include"].map({True: "black", False: "red"})

plot_embedding(
    ax=axs[1, 0],
    embedding=principal_components,
    x='PC1',
    y='PC2',
    z='CellCount',
    thresholds=COUNTS_THRESHOLDS,
    title='PCA Norm. Counts\nColored by Cell Count',
    include_colorbar=True
)
plot_embedding(
    ax=axs[1, 1],
    embedding=principal_components,
    x='PC1',
    y='PC2',
    z='hue',
    thresholds=COUNTS_THRESHOLDS,
    title='PCA Norm. Counts\nColored by MAD outlier',
    include_colorbar=False
)
plot_legend(
    ax=axs[1, 2],
    df=principal_components,
    x='PC1',
    y='PC2',
    thresholds=COUNTS_THRESHOLDS,
)

include_principal_components = filter_on_mad(df=principal_components, x="PC1", y="PC2", threshold=args.counts_threshold)
include_principal_components = include_principal_components.loc[:, ["include"]].copy()
include_principal_components.columns = ["counts"]

smf = smf.merge(include_principal_components, left_on=sample_aggregate, right_index=True, how="left")

################################################################################

print("Processing full kinship")
full_kinship_pca = pd.read_csv(args.full_kinship, sep="\t", header=0, index_col=0, low_memory=False)
full_kinship_pca["include"] = "red"
full_kinship_pca.loc[[sample for sample in smf[individual_aggregate] if sample in full_kinship_pca.index], "include"] = "black"

plot_embedding(
    ax=axs[2, 0],
    embedding=full_kinship_pca,
    x='PC1',
    y='PC2',
    z='include',
    alpha=0.6,
    thresholds=KINSHIP_THRESHOLDS,
    title='Kinship\nAll samples',
    include_colorbar=False
)
x_data = full_kinship_pca.loc[full_kinship_pca["include"] == "black", "PC1"].to_numpy()
y_data = full_kinship_pca.loc[full_kinship_pca["include"] == "black", "PC2"].to_numpy()
plot_scatterplot(ax=axs[2, 0], x_data=x_data, y_data=y_data)
del (x_data, y_data)

print("Processing kinship")
kinship_pca = pd.read_csv(args.kinship, sep="\t", header=0, index_col=0, low_memory=False)
kinship_pca = filter_on_mad(df=kinship_pca, x="PC1", y="PC2", threshold=args.kinship_threshold)
kinship_pca["hue"] = kinship_pca["include"].map({True: "black", False: "red"})

plot_embedding(
    ax=axs[2, 1],
    embedding=kinship_pca,
    x='PC1',
    y='PC2',
    z='hue',
    thresholds=KINSHIP_THRESHOLDS,
    title='Kinship\nIncluded samples',
    include_colorbar=False
)
plot_legend(
    ax=axs[2, 2],
    df=kinship_pca,
    x='PC1',
    y='PC2',
    thresholds=KINSHIP_THRESHOLDS,
)

include_kinship_pca = filter_on_mad(df=kinship_pca, x="PC1", y="PC2", threshold=args.kinship_threshold)
include_kinship_pca = include_kinship_pca.loc[:, ["include"]].copy()
include_kinship_pca.columns = ["kinship"]

smf = smf.merge(include_kinship_pca, left_on=individual_aggregate, right_index=True, how="left")

################################################################################

filter_columns = [col for col in smf.columns if col not in [individual_aggregate, sample_aggregate]]
smf["all"] = smf[filter_columns].sum(axis=1) == len(filter_columns)
n_preqc = smf.shape[0]
n_postqc = smf["all"].sum()

fig.suptitle("{} - {} - {}\n Pre-filter N = {:,}    Post-filter N = {:,}".format(args.ancestry, args.cell_level, args.cell_type, n_preqc, n_postqc))
fig.tight_layout()
plt.savefig(args.plot_out + 'sample_qc.png', bbox_inches="tight")
print("")

smf.to_csv(args.data_out + 'sample_qc.txt', sep="\t", header=True, index=False)

################################################################################

print("Filtering samples")
for column in filter_columns:
    print("    {:,} samples passed {} thresholds".format(np.sum(smf[column]), column))
smf = smf.loc[~smf["all"], [individual_aggregate, sample_aggregate]]
print("    {:,} samples failed QC thresholds".format(smf.shape[0]))
smf.to_csv(args.data_out + 'exclude_smf.txt', sep="\t", header=False, index=False)

print("Saving filter settings")
with gzopen(args.data_out + 'threshold_selection.txt', mode='w') as f:
    f.write("Threshold\tValue\tPASS\n")
    f.write("donor_cell_threshold\t{}\tFalse\n".format(args.donor_cell_threshold))
    f.write("counts_threshold\t{}\tFalse\n".format(args.counts_threshold))
    f.write("kinship_threshold\t{}\tFalse\n".format(args.kinship_threshold))
f.close()
