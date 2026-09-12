# Copyright (C) <2019> University of Bern - Interfaculty Bioinformatics Unit
#
# This program is heavily based on Marco Galardini's roary_plots.
# https://github.com/sanger-pathogens/Roary/tree/master/contrib/roary_plots
#
# Databiomics publication-ready extensions preserve the original workflow while
# adding multi-format export and additional comparative-genomics visualisations.

__author__ = "Thomas Roder; Databiomics extensions"
__version__ = "0.3.0"

import json
import math
import os
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import Phylo

from .utils import load_og, load_hog


ALLOWED_FORMATS = ("png", "tiff", "pdf", "svg")
# Okabe-Ito inspired, colour-vision-deficiency friendly palette.
COLORS = {
    "blue": "#0072B2",
    "sky": "#56B4E9",
    "green": "#009E73",
    "orange": "#E69F00",
    "vermillion": "#D55E00",
    "purple": "#CC79A7",
    "grey": "#7A7A7A",
    "light_grey": "#E6E6E6",
    "dark": "#202124",
}


def import_tree(path_to_newick):
    path_to_newick = os.path.abspath(os.path.expanduser(str(path_to_newick)))
    if not os.path.isfile(path_to_newick):
        raise FileNotFoundError(f'Newick tree does not exist: "{path_to_newick}"')
    return Phylo.read(path_to_newick, "newick")


def import_roary_table(path_to_table, skipped_columns=14):
    """Backwards-compatible helper for Roary presence/absence CSV files."""
    if not os.path.isfile(path_to_table):
        raise FileNotFoundError(path_to_table)
    pandas_table = pd.read_csv(path_to_table, sep=",", low_memory=False)
    pandas_table.set_index("Gene", inplace=True)
    pandas_table.drop(list(pandas_table.columns[:skipped_columns - 1]), axis=1, inplace=True)
    pandas_table.replace(".{2,100}", 1, regex=True, inplace=True)
    pandas_table.replace(np.nan, 0, regex=True, inplace=True)
    return pandas_table.transpose()


def _parse_formats(format="svg", formats=None):
    values = formats if formats is not None else format
    if isinstance(values, str):
        values = [part.strip().lower().lstrip(".") for part in values.split(",") if part.strip()]
    else:
        values = [str(part).strip().lower().lstrip(".") for part in values]

    if not values:
        values = ["svg"]

    invalid = sorted(set(values) - set(ALLOWED_FORMATS))
    if invalid:
        raise ValueError(f"Unsupported figure format(s): {', '.join(invalid)}. Allowed: {', '.join(ALLOWED_FORMATS)}")

    # Preserve caller order while removing duplicates.
    return list(dict.fromkeys(values))


def _publication_style():
    sns.set_theme(context="paper", style="white", font_scale=1.05)
    plt.rcParams.update({
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "axes.edgecolor": COLORS["dark"],
        "axes.labelcolor": COLORS["dark"],
        "text.color": COLORS["dark"],
        "xtick.color": COLORS["dark"],
        "ytick.color": COLORS["dark"],
        "axes.titleweight": "bold",
        "axes.titlesize": 12,
        "axes.labelsize": 10,
        "savefig.facecolor": "white",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })


def _save_figure(fig, out, stem, formats, dpi):
    generated = []
    for fmt in formats:
        path = Path(out) / f"{stem}.{fmt}"
        kwargs = {"bbox_inches": "tight", "pad_inches": 0.08}
        if fmt in {"png", "tiff"}:
            kwargs["dpi"] = dpi
        fig.savefig(path, **kwargs)
        generated.append(str(path))
    plt.close(fig)
    return generated


def _presence_categories(og_count, n_genomes):
    core_threshold = max(1, math.ceil(n_genomes * 0.99))
    soft_threshold = max(1, math.ceil(n_genomes * 0.95))
    shell_threshold = max(1, math.ceil(n_genomes * 0.15))

    core = int(((og_count >= core_threshold) & (og_count <= n_genomes)).sum())
    softcore = int(((og_count >= soft_threshold) & (og_count < core_threshold)).sum())
    shell = int(((og_count >= shell_threshold) & (og_count < soft_threshold)).sum())
    cloud = int((og_count < shell_threshold).sum())

    return {
        "core": core,
        "softcore": softcore,
        "shell": shell,
        "cloud": cloud,
        "thresholds": {
            "core_min_genomes": core_threshold,
            "softcore_min_genomes": soft_threshold,
            "shell_min_genomes": shell_threshold,
        },
    }


def _tree_order(tree, columns):
    tips = [terminal.name for terminal in tree.get_terminals()]
    missing = [name for name in tips if name not in columns]
    if missing:
        preview = ", ".join(str(x) for x in missing[:8])
        suffix = "..." if len(missing) > 8 else ""
        raise ValueError(
            "Species-tree tip labels must match Orthogroups.tsv column names. "
            f"Missing {len(missing)} tip(s) from the table: {preview}{suffix}"
        )
    return tips


def _plot_frequency(og_count, n_genomes):
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    bins = np.arange(0.5, n_genomes + 1.5, 1)
    ax.hist(
        og_count.to_numpy(dtype=float),
        bins=bins,
        color=COLORS["blue"],
        edgecolor="white",
        linewidth=0.5,
        alpha=0.9,
    )
    ax.set_xlabel("Number of genomes containing the orthogroup")
    ax.set_ylabel("Number of orthogroups")
    ax.set_title("Orthogroup frequency across genomes")
    ax.set_xlim(0.5, n_genomes + 0.5)
    if n_genomes <= 30:
        ax.set_xticks(np.arange(1, n_genomes + 1))
    sns.despine(ax=ax)
    fig.tight_layout()
    return fig


def _plot_matrix(tree, matrix, no_labels=False):
    n_orthogroups, n_genomes = matrix.shape
    tips = _tree_order(tree, matrix.columns)
    ordered = matrix.loc[:, tips]

    # Width increases mildly for large datasets while keeping manuscript-friendly proportions.
    matrix_width = min(14.0, max(7.0, 5.0 + n_orthogroups / 800.0))
    tree_width = min(6.5, max(3.5, 3.8 + n_genomes / 80.0))
    fig_height = min(18.0, max(6.0, 3.5 + n_genomes * 0.22))

    fig = plt.figure(figsize=(tree_width + matrix_width, fig_height))
    grid = fig.add_gridspec(1, 2, width_ratios=[tree_width, matrix_width], wspace=0.01)
    ax_tree = fig.add_subplot(grid[0, 0])
    ax_matrix = fig.add_subplot(grid[0, 1])

    label_size = max(5.0, min(9.0, 11.0 - n_genomes * 0.06))
    Phylo.draw(
        tree,
        axes=ax_tree,
        show_confidence=False,
        label_func=(lambda _: None) if no_labels else (lambda node: str(node.name or "")),
        do_show=False,
    )
    ax_tree.set_title(f"Species tree ({n_genomes} genomes)")
    ax_tree.set_xlabel("")
    ax_tree.set_ylabel("")
    ax_tree.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
    for text in ax_tree.texts:
        text.set_fontsize(label_size)

    data = ordered.T.to_numpy(dtype=np.uint8)
    cmap = matplotlib.colors.ListedColormap(["#F2F4F7", COLORS["blue"]])
    ax_matrix.imshow(data, cmap=cmap, vmin=0, vmax=1, aspect="auto", interpolation="nearest")
    ax_matrix.set_title(f"Orthogroup presence/absence ({n_orthogroups:,} clusters)")
    ax_matrix.set_xlabel("Orthogroups (ordered by prevalence)")
    ax_matrix.set_ylabel("")
    ax_matrix.set_yticks([])
    ax_matrix.set_xticks([])
    for spine in ax_matrix.spines.values():
        spine.set_visible(False)

    fig.subplots_adjust(left=0.04, right=0.995, top=0.93, bottom=0.05)
    return fig


def _plot_composition(categories, n_orthogroups):
    labels = ["Core", "Soft-core", "Shell", "Cloud"]
    values = [categories["core"], categories["softcore"], categories["shell"], categories["cloud"]]
    palette = [COLORS["blue"], COLORS["green"], COLORS["orange"], COLORS["grey"]]

    fig, ax = plt.subplots(figsize=(7.2, 6.2))
    wedges, _ = ax.pie(
        values,
        startangle=90,
        counterclock=False,
        colors=palette,
        wedgeprops={"width": 0.42, "edgecolor": "white", "linewidth": 1.2},
    )
    ax.text(0, 0.06, f"{n_orthogroups:,}", ha="center", va="center", fontsize=18, fontweight="bold")
    ax.text(0, -0.10, "orthogroups", ha="center", va="center", fontsize=9)

    legend_labels = []
    for label, value in zip(labels, values):
        pct = 100.0 * value / n_orthogroups if n_orthogroups else 0.0
        legend_labels.append(f"{label}: {value:,} ({pct:.1f}%)")
    ax.legend(wedges, legend_labels, loc="center left", bbox_to_anchor=(0.96, 0.5), frameon=False)
    ax.set_title("Pangenome composition")
    fig.tight_layout()
    return fig


def _accumulation_statistics(matrix, permutations=100, seed=42):
    arr = matrix.to_numpy(dtype=np.uint8)
    n_orthogroups, n_genomes = arr.shape
    if n_genomes == 0:
        raise ValueError("No genomes were found in the orthogroup table.")

    permutations = max(1, min(int(permutations), 500))
    if n_genomes == 1:
        pan = np.array([[int(arr[:, 0].sum())]], dtype=float)
        core = pan.copy()
    else:
        rng = np.random.default_rng(seed)
        pan = np.zeros((permutations, n_genomes), dtype=float)
        core = np.zeros((permutations, n_genomes), dtype=float)
        for idx in range(permutations):
            order = rng.permutation(n_genomes)
            permuted = arr[:, order]
            pan[idx, :] = np.maximum.accumulate(permuted, axis=1).sum(axis=0)
            core[idx, :] = np.minimum.accumulate(permuted, axis=1).sum(axis=0)

    def stats(values):
        return {
            "mean": values.mean(axis=0),
            "low": np.percentile(values, 2.5, axis=0),
            "high": np.percentile(values, 97.5, axis=0),
        }

    return stats(pan), stats(core)


def _plot_accumulation(matrix, permutations=100):
    pan, core = _accumulation_statistics(matrix, permutations=permutations)
    x = np.arange(1, matrix.shape[1] + 1)

    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    ax.plot(x, pan["mean"], color=COLORS["blue"], linewidth=2.2, label="Pan-genome")
    ax.fill_between(x, pan["low"], pan["high"], color=COLORS["blue"], alpha=0.15, linewidth=0)
    ax.plot(x, core["mean"], color=COLORS["vermillion"], linewidth=2.2, label="Core genome")
    ax.fill_between(x, core["low"], core["high"], color=COLORS["vermillion"], alpha=0.15, linewidth=0)
    ax.set_xlabel("Number of genomes sampled")
    ax.set_ylabel("Number of orthogroups")
    ax.set_title(f"Pan/core-genome accumulation ({max(1, int(permutations))} random orders)")
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.18)
    sns.despine(ax=ax)
    fig.tight_layout()
    return fig


def _top_patterns(matrix, max_patterns=20):
    max_patterns = max(4, min(int(max_patterns), 50))
    columns = list(matrix.columns)
    counts = Counter(tuple(row) for row in matrix.to_numpy(dtype=np.uint8))
    top = counts.most_common(max_patterns)
    patterns = np.asarray([pattern for pattern, _ in top], dtype=np.uint8)
    sizes = np.asarray([count for _, count in top], dtype=int)
    return columns, patterns, sizes


def _plot_patterns(matrix, max_patterns=20):
    columns, patterns, sizes = _top_patterns(matrix, max_patterns=max_patterns)
    n_patterns = len(sizes)
    fig_height = min(12.0, max(5.0, 2.8 + 0.34 * n_patterns))
    fig = plt.figure(figsize=(12.5, fig_height))
    grid = fig.add_gridspec(1, 2, width_ratios=[3.2, 6.8], wspace=0.08)
    ax_bar = fig.add_subplot(grid[0, 0])
    ax_mat = fig.add_subplot(grid[0, 1])

    y = np.arange(n_patterns)
    ax_bar.barh(y, sizes, color=COLORS["blue"], alpha=0.9)
    ax_bar.invert_yaxis()
    ax_bar.set_xlabel("Orthogroups")
    ax_bar.set_ylabel("Presence/absence pattern")
    ax_bar.set_yticks(y)
    ax_bar.set_yticklabels([f"Pattern {i + 1}" for i in y], fontsize=8)
    ax_bar.grid(axis="x", alpha=0.16)
    sns.despine(ax=ax_bar)

    cmap = matplotlib.colors.ListedColormap(["#F2F4F7", COLORS["blue"]])
    ax_mat.imshow(patterns, cmap=cmap, vmin=0, vmax=1, aspect="auto", interpolation="nearest")
    ax_mat.set_yticks(y)
    ax_mat.set_yticklabels([])
    if len(columns) <= 40:
        ax_mat.set_xticks(np.arange(len(columns)))
        ax_mat.set_xticklabels(columns, rotation=90, fontsize=max(5.0, 8.0 - len(columns) * 0.04))
    else:
        ax_mat.set_xticks([])
        ax_mat.set_xlabel(f"{len(columns)} genomes (columns follow Orthogroups.tsv order)")
    ax_mat.set_title(f"Top {n_patterns} orthogroup occupancy patterns")
    for spine in ax_mat.spines.values():
        spine.set_visible(False)

    fig.subplots_adjust(left=0.08, right=0.995, top=0.92, bottom=0.16, wspace=0.08)
    return fig


def create_plots(
    tree,
    orthogroups_tsv,
    out,
    format="svg",
    no_labels=False,
    hog=False,
    formats=None,
    dpi=600,
    permutations=100,
    max_patterns=20,
):
    """
    Create publication-ready comparative-genomics figures from OrthoFinder OG/HOG outputs.

    Backwards compatibility:
      * ``format='svg'`` still works exactly as the legacy CLI expects.
      * Set ``formats='svg,pdf,png'`` to export multiple formats in one run.

    Parameters
    ----------
    tree
        Newick tree object or path to a Newick species tree.
    orthogroups_tsv
        Path to Orthogroups.tsv/N0.tsv or a boolean DataFrame
        (rows=orthogroups, columns=genomes).
    out
        Output directory.
    format
        Legacy single output format.
    no_labels
        Hide species labels on the phylogenetic tree.
    hog
        If True, read a hierarchical orthogroup N0.tsv file.
    formats
        Optional comma-separated or iterable list of png,tiff,pdf,svg.
    dpi
        Raster export resolution. Defaults to 600 dpi.
    permutations
        Random genome orders used for pan/core accumulation confidence bands.
    max_patterns
        Maximum exact presence/absence patterns shown in the pattern panel.
    """
    out = os.path.abspath(os.path.expanduser(str(out)))
    os.makedirs(out, exist_ok=True)
    export_formats = _parse_formats(format=format, formats=formats)
    dpi = max(150, min(int(dpi), 1200))

    if isinstance(tree, (str, os.PathLike)):
        tree = import_tree(tree)
    if not isinstance(tree, Phylo.Newick.Tree):
        raise TypeError("tree must be a Bio.Phylo Newick Tree or a path to a Newick file")

    if isinstance(orthogroups_tsv, (str, os.PathLike)):
        if hog:
            orthogroups_tsv = load_hog(str(orthogroups_tsv), result_type="boolean")
        else:
            orthogroups_tsv = load_og(str(orthogroups_tsv), result_type="boolean")
    if not isinstance(orthogroups_tsv, pd.DataFrame):
        raise TypeError("orthogroups_tsv must be a path or pandas DataFrame")

    matrix = orthogroups_tsv.astype(bool)
    if matrix.empty or matrix.shape[1] == 0:
        raise ValueError("Orthogroup table is empty.")

    _publication_style()
    n_orthogroups, n_genomes = matrix.shape
    og_count = matrix.sum(axis=1)
    prevalence_order = og_count.sort_values(ascending=False).index
    matrix_sorted = matrix.loc[prevalence_order]
    categories = _presence_categories(og_count, n_genomes)

    outputs = []
    outputs.extend(_save_figure(_plot_frequency(og_count, n_genomes), out, "pangenome_frequency", export_formats, dpi))
    outputs.extend(_save_figure(_plot_matrix(tree, matrix_sorted, no_labels=no_labels), out, "pangenome_matrix", export_formats, dpi))
    # Keep the historic filename while upgrading the visual to a publication-ready donut.
    outputs.extend(_save_figure(_plot_composition(categories, n_orthogroups), out, "pangenome_pie", export_formats, dpi))
    outputs.extend(_save_figure(_plot_accumulation(matrix, permutations=permutations), out, "pangenome_accumulation", export_formats, dpi))
    outputs.extend(_save_figure(_plot_patterns(matrix, max_patterns=max_patterns), out, "orthogroup_patterns", export_formats, dpi))

    summary = {
        "schema_version": 1,
        "n_genomes": int(n_genomes),
        "n_orthogroups": int(n_orthogroups),
        "core": categories["core"],
        "softcore": categories["softcore"],
        "shell": categories["shell"],
        "cloud": categories["cloud"],
        "thresholds": categories["thresholds"],
        "formats": export_formats,
        "raster_dpi": dpi,
        "accumulation_permutations": max(1, min(int(permutations), 500)),
        "outputs": [os.path.basename(path) for path in outputs],
    }
    summary_path = Path(out) / "orthofinder_plot_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    outputs.append(str(summary_path))

    print("OrthoFinder plot summary:", json.dumps(summary, sort_keys=True))
    return summary


def main():
    import fire
    fire.Fire(create_plots)


if __name__ == "__main__":
    main()
