"""
Synchronicity calculation and plotting functions.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Dict, List, Optional, Callable

from utils import (
    average_pairwise_cosine_similarity,
    average_pairwise_euclidean_similarity,
    bootstrap_confidence_interval,
)


def compute_synchronicity(
    adatas: Dict[str, "AnnData"],
    adatas_acaA: Dict[str, "AnnData"],
    adatas_acaA_pkaC: Dict[str, "AnnData"],
    time_points: List[str],
    embedding_key: str = "X_uce",
    similarity_func: Callable = average_pairwise_euclidean_similarity,
    ax4_markers: List[str] = ["act15GFP", "mCerulean", "mCherry"],
    acaA_markers: List[str] = ["act15GFP", "mCerulean", "mCherry", "mNeonG"],
    acaA_pkaC_markers: List[str] = ["r1", "r2", "r3"],
    late_timepoints: List[str] = ["12hr", "16hr", "20hr"],
) -> Dict:
    """
    Compute synchronicity metrics for all strains, timepoints, and subsets.

    Parameters
    ----------
    adatas : dict
        Dictionary mapping timepoint -> AnnData for AX4
    adatas_acaA : dict
        Dictionary mapping timepoint -> AnnData for acaA
    adatas_acaA_pkaC : dict
        Dictionary mapping timepoint -> AnnData for acaA_pkaC
    time_points : list
        List of timepoint strings, e.g. ["00hr", "04hr", ...]
    embedding_key : str
        Key in adata.obsm to use for cell vectors (e.g. "X_pca", "X_uce")
    similarity_func : callable
        Function to compute similarity from cell vectors (e.g. average_pairwise_cosine_similarity)
    ax4_markers : list
        Marker names for AX4 replicates
    acaA_markers : list
        Marker names for acaA replicates
    acaA_pkaC_markers : list
        Marker names for acaA_pkaC replicates
    late_timepoints : list
        Timepoints where prespore/prestalk cells exist

    Returns
    -------
    dict
        Dictionary with keys:
        - 'mean_distances': synchronicity values
        - 'cell_counts': cell counts
        - 'median_genes': median genes per cell
        - 'total_reads': total reads
    """

    # Initialize result dictionaries
    mean_distances = {
        "AX4": {"all": {}},
        "AX4_psp": {"all": {}},
        "AX4_pst": {"all": {}},
        "acaA": {"all": {}},
        "acaA_pkaC": {"all": {}},
        "acaA_pkaC_psp": {"all": {}},
        "acaA_pkaC_pst": {"all": {}},
    }

    cell_counts = {
        "AX4": {"all": {}},
        "AX4_psp": {"all": {}},
        "AX4_pst": {"all": {}},
        "acaA": {"all": {}},
        "acaA_pkaC": {"all": {}},
        "acaA_pkaC_psp": {"all": {}},
        "acaA_pkaC_pst": {"all": {}},
    }

    median_genes = {
        "AX4": {"all": {}},
        "AX4_psp": {"all": {}},
        "AX4_pst": {"all": {}},
        "acaA": {"all": {}},
        "acaA_pkaC": {"all": {}},
        "acaA_pkaC_psp": {"all": {}},
        "acaA_pkaC_pst": {"all": {}},
    }

    total_reads = {
        "AX4": {"all": {}},
        "AX4_psp": {"all": {}},
        "AX4_pst": {"all": {}},
        "acaA": {"all": {}},
        "acaA_pkaC": {"all": {}},
        "acaA_pkaC_psp": {"all": {}},
        "acaA_pkaC_pst": {"all": {}},
    }

    # Add marker keys
    for marker in ax4_markers:
        mean_distances["AX4"][marker] = {}
        mean_distances["AX4_psp"][marker] = {}
        mean_distances["AX4_pst"][marker] = {}
        cell_counts["AX4"][marker] = {}
        cell_counts["AX4_psp"][marker] = {}
        cell_counts["AX4_pst"][marker] = {}
        median_genes["AX4"][marker] = {}
        median_genes["AX4_psp"][marker] = {}
        median_genes["AX4_pst"][marker] = {}
        total_reads["AX4"][marker] = {}
        total_reads["AX4_psp"][marker] = {}
        total_reads["AX4_pst"][marker] = {}

    for marker in acaA_markers:
        mean_distances["acaA"][marker] = {}
        cell_counts["acaA"][marker] = {}
        median_genes["acaA"][marker] = {}
        total_reads["acaA"][marker] = {}

    for marker in acaA_pkaC_markers:
        mean_distances["acaA_pkaC"][marker] = {}
        mean_distances["acaA_pkaC_psp"][marker] = {}
        mean_distances["acaA_pkaC_pst"][marker] = {}
        cell_counts["acaA_pkaC"][marker] = {}
        cell_counts["acaA_pkaC_psp"][marker] = {}
        cell_counts["acaA_pkaC_pst"][marker] = {}
        median_genes["acaA_pkaC"][marker] = {}
        median_genes["acaA_pkaC_psp"][marker] = {}
        median_genes["acaA_pkaC_pst"][marker] = {}
        total_reads["acaA_pkaC"][marker] = {}
        total_reads["acaA_pkaC_psp"][marker] = {}
        total_reads["acaA_pkaC_pst"][marker] = {}

    # Helper function to compute metrics for a subset
    def compute_metrics(adata_subset, strain_key, subset_key, tp):
        if embedding_key not in adata_subset.obsm:
            return

        cell_vectors = adata_subset.obsm[embedding_key]
        if len(cell_vectors) < 2:
            return

        # Compute n_genes and total_reads if not present
        X = adata_subset.X
        if hasattr(X, "toarray"):
            X = X.toarray()
        n_genes = (X > 0).sum(axis=1)
        total_reads_vals = X.sum(axis=1)

        mean_distances[strain_key][subset_key][tp] = similarity_func(cell_vectors)
        cell_counts[strain_key][subset_key][tp] = len(cell_vectors)
        median_genes[strain_key][subset_key][tp] = np.median(n_genes)
        total_reads[strain_key][subset_key][tp] = np.sum(total_reads_vals)

    # Process AX4
    for tp in time_points:
        adata = adatas[tp]

        # All cells
        compute_metrics(adata, "AX4", "all", tp)

        # Per marker
        for marker in ax4_markers:
            subset = adata[adata.obs["marker"] == marker]
            if len(subset) > 0:
                compute_metrics(subset, "AX4", marker, tp)

        # Prespore/prestalk (only for late timepoints)
        if tp in late_timepoints and "cell_type" in adata.obs.columns:
            for cell_type, ct_key in [("psp", "AX4_psp"), ("pst", "AX4_pst")]:
                subset = adata[adata.obs["cell_type"].str.contains(cell_type, na=False)]
                if len(subset) > 0:
                    compute_metrics(subset, ct_key, "all", tp)

                    # Per marker within cell type
                    for marker in ax4_markers:
                        subset2 = subset[subset.obs["marker"] == marker]
                        if len(subset2) > 0:
                            compute_metrics(subset2, ct_key, marker, tp)

    # Process acaA
    for tp in time_points:
        adata = adatas_acaA[tp]

        # All cells
        compute_metrics(adata, "acaA", "all", tp)

        # Per marker
        for marker in acaA_markers:
            subset = adata[adata.obs["marker"] == marker]
            if len(subset) > 0:
                compute_metrics(subset, "acaA", marker, tp)

    # Process acaA_pkaC
    for tp in time_points:
        adata = adatas_acaA_pkaC[tp]

        # All cells
        compute_metrics(adata, "acaA_pkaC", "all", tp)

        # Per marker
        for marker in acaA_pkaC_markers:
            subset = adata[adata.obs["marker"] == marker]
            if len(subset) > 0:
                compute_metrics(subset, "acaA_pkaC", marker, tp)

        # Prespore/prestalk (only for late timepoints)
        if tp in late_timepoints and "cell_type" in adata.obs.columns:
            for cell_type, ct_key in [
                ("psp", "acaA_pkaC_psp"),
                ("pst", "acaA_pkaC_pst"),
            ]:
                subset = adata[adata.obs["cell_type"].str.contains(cell_type, na=False)]
                if len(subset) > 0:
                    compute_metrics(subset, ct_key, "all", tp)

                    # Per marker within cell type
                    for marker in acaA_pkaC_markers:
                        subset2 = subset[subset.obs["marker"] == marker]
                        if len(subset2) > 0:
                            compute_metrics(subset2, ct_key, marker, tp)

    return {
        "mean_distances": mean_distances,
        "cell_counts": cell_counts,
        "median_genes": median_genes,
        "total_reads": total_reads,
    }


def compute_confidence_intervals(
    mean_distances: Dict,
    time_points: List[str],
    ax4_markers: List[str] = ["act15GFP", "mCerulean", "mCherry"],
    acaA_markers: List[str] = ["act15GFP", "mCerulean", "mCherry", "mNeonG"],
    acaA_pkaC_markers: List[str] = ["r1", "r2", "r3"],
    late_timepoints: List[str] = ["12hr", "16hr", "20hr"],
    num_bootstrap: int = 10000,
    ci: int = 95,
    random_state: Optional[int] = None,
) -> Dict:
    """
    Compute bootstrap confidence intervals for synchronicity across markers/replicates.

    Parameters
    ----------
    mean_distances : dict
        Output from compute_synchronicity['mean_distances']
    time_points : list
        List of timepoint strings, e.g. ["00hr", "04hr", ...]
    ax4_markers : list
        Marker names for AX4 replicates
    acaA_markers : list
        Marker names for acaA replicates
    acaA_pkaC_markers : list
        Marker names for acaA_pkaC replicates
    late_timepoints : list
        Timepoints where prespore/prestalk cells exist
    num_bootstrap : int
        Number of bootstrap samples (default 10000)
    ci : int
        Confidence interval percentage (default 95)
    random_state : int, optional
        Random seed for reproducible bootstrap results

    Returns
    -------
    dict
        Dictionary with keys for each strain/cell type, containing:
        {timepoint: [mean, lower, upper, bootstrap_means]}
    """
    if random_state is not None:
        np.random.seed(random_state)
    CI_dict = {
        "AX4": {},
        "acaA": {},
        "acaA_pkaC": {},
        "AX4_psp": {},
        "AX4_pst": {},
        "acaA_pkaC_psp": {},
        "acaA_pkaC_pst": {},
    }

    # AX4 confidence intervals
    for tp in time_points:
        # All cells - bootstrap across markers
        bootstrap_means = []
        for marker in ax4_markers:
            if marker in mean_distances["AX4"] and tp in mean_distances["AX4"][marker]:
                bootstrap_means.append(mean_distances["AX4"][marker][tp])
        if bootstrap_means:
            CI_dict["AX4"][tp] = bootstrap_confidence_interval(
                bootstrap_means, num_bootstrap=num_bootstrap, ci=ci
            )

        # Prespore/prestalk (only for late timepoints)
        if tp in late_timepoints:
            for cell_type in ["psp", "pst"]:
                strain_key = f"AX4_{cell_type}"
                bootstrap_means = []
                for marker in ax4_markers:
                    if (
                        marker in mean_distances[strain_key]
                        and tp in mean_distances[strain_key][marker]
                    ):
                        bootstrap_means.append(mean_distances[strain_key][marker][tp])
                if bootstrap_means:
                    CI_dict[strain_key][tp] = bootstrap_confidence_interval(
                        bootstrap_means, num_bootstrap=num_bootstrap, ci=ci
                    )

    # acaA confidence intervals
    for tp in time_points:
        bootstrap_means = []
        for marker in acaA_markers:
            if (
                marker in mean_distances["acaA"]
                and tp in mean_distances["acaA"][marker]
            ):
                bootstrap_means.append(mean_distances["acaA"][marker][tp])
        if bootstrap_means:
            CI_dict["acaA"][tp] = bootstrap_confidence_interval(
                bootstrap_means, num_bootstrap=num_bootstrap, ci=ci
            )

    # acaA_pkaC confidence intervals
    for tp in time_points:
        # All cells - bootstrap across markers
        bootstrap_means = []
        for marker in acaA_pkaC_markers:
            if (
                marker in mean_distances["acaA_pkaC"]
                and tp in mean_distances["acaA_pkaC"][marker]
            ):
                bootstrap_means.append(mean_distances["acaA_pkaC"][marker][tp])
        if bootstrap_means:
            CI_dict["acaA_pkaC"][tp] = bootstrap_confidence_interval(
                bootstrap_means, num_bootstrap=num_bootstrap, ci=ci
            )

        # Prespore/prestalk (only for late timepoints)
        if tp in late_timepoints:
            for cell_type in ["psp", "pst"]:
                strain_key = f"acaA_pkaC_{cell_type}"
                bootstrap_means = []
                for marker in acaA_pkaC_markers:
                    if (
                        marker in mean_distances[strain_key]
                        and tp in mean_distances[strain_key][marker]
                    ):
                        bootstrap_means.append(mean_distances[strain_key][marker][tp])
                if bootstrap_means:
                    CI_dict[strain_key][tp] = bootstrap_confidence_interval(
                        bootstrap_means, num_bootstrap=num_bootstrap, ci=ci
                    )

    return CI_dict


def plot_synchronicity_row(
    mean_distances: Dict,
    ax4_markers: List[str] = ["act15GFP", "mCerulean", "mCherry", "mNeonG"],
    acaA_markers: List[str] = ["act15GFP", "mCerulean", "mCherry", "mNeonG"],
    acaA_pkaC_markers: List[str] = ["r1", "r2", "r3"],
    suptitle: Optional[str] = None,
    save_path: Optional[str] = None,
    figsize: tuple = (24.462, 6.257),
    share_y: bool = False,
    CI_dict: Optional[Dict] = None,
    capsize: int = 5,
    capthick: int = 2,
):
    """
    Plot all strain synchronicity plots side by side in a single row.

    Parameters
    ----------
    mean_distances : dict
        Output from compute_synchronicity['mean_distances']
    ax4_markers, acaA_markers, acaA_pkaC_markers : list
        Marker lists for each strain
    suptitle : str, optional
        Super title for the entire figure
    save_path : str, optional
        Path to save figure
    figsize : tuple
        Figure size (width, height)
    share_y : bool
        Whether to share y-axis across all subplots
    CI_dict : dict, optional
        Confidence interval dictionary from compute_confidence_intervals.
        If provided, error bars will be drawn for "all" cell lines.
    capsize : int
        Size of the error bar caps (default 5)
    capthick : int
        Thickness of the error bar caps (default 2)

    Returns
    -------
    fig, axes
    """
    fig, axes = plt.subplots(1, 4, figsize=figsize, sharey=share_y)

    # Helper to convert timepoint labels (e.g., "00hr" -> "0", "04hr" -> "4")
    def format_timepoint(tp):
        return str(int(tp.replace("hr", "")))

    marker_colors = {
        "act15GFP": "#3b76b0",
        "mCerulean": "orange",
        "mCherry": "green",
        "mNeonG": "#8d69b8",
        "r1": "#3b76b0",
        "r2": "orange",
        "r3": "green",
    }

    marker_labels = {
        "act15GFP": "GFP",
        "mCerulean": "mCerulean",
        "mCherry": "mCherry",
        "mNeonG": "mNeonG",
        "r1": "replicate 1",
        "r2": "replicate 2",
        "r3": "replicate 3",
    }

    # Helper function to add error bars
    def add_error_bars(ax, strain_key, time_points_list, color, linestyle="-"):
        if CI_dict is None or strain_key not in CI_dict:
            return
        # Filter to only timepoints that have CI data
        valid_tps = [tp for tp in time_points_list if tp in CI_dict[strain_key]]
        if not valid_tps:
            return
        yerr_lower = [
            CI_dict[strain_key][tp][0] - CI_dict[strain_key][tp][1] for tp in valid_tps
        ]
        yerr_upper = [
            CI_dict[strain_key][tp][2] - CI_dict[strain_key][tp][0] for tp in valid_tps
        ]
        yvals = [mean_distances[strain_key]["all"][tp] for tp in valid_tps]

        ax.errorbar(
            valid_tps,
            yvals,
            yerr=[yerr_lower, yerr_upper],
            fmt="none",
            color=color,
            capsize=capsize,
            capthick=capthick,
            linestyle=linestyle,
        )

    # Panel 1: All strains combined
    ax = axes[0]
    if mean_distances["AX4"]["all"]:
        tps = list(mean_distances["AX4"]["all"].keys())
        ax.plot(
            tps,
            list(mean_distances["AX4"]["all"].values()),
            label="AX4",
            color="black",
            marker="o",
        )
        add_error_bars(ax, "AX4", tps, "black")
    if mean_distances["AX4_psp"]["all"]:
        tps = list(mean_distances["AX4_psp"]["all"].keys())
        ax.plot(
            tps,
            list(mean_distances["AX4_psp"]["all"].values()),
            label="AX4 prespore",
            color="deepskyblue",
            marker="o",
        )
        add_error_bars(ax, "AX4_psp", tps, "deepskyblue")
    if mean_distances["AX4_pst"]["all"]:
        tps = list(mean_distances["AX4_pst"]["all"].keys())
        ax.plot(
            tps,
            list(mean_distances["AX4_pst"]["all"].values()),
            label="AX4 prestalk",
            color="red",
            marker="o",
        )
        add_error_bars(ax, "AX4_pst", tps, "red")
    if mean_distances["acaA"]["all"]:
        tps = list(mean_distances["acaA"]["all"].keys())
        ax.plot(
            tps,
            list(mean_distances["acaA"]["all"].values()),
            label="acaA−",
            color="black",
            marker="o",
            linestyle="dashed",
        )
        add_error_bars(ax, "acaA", tps, "black")
    if mean_distances["acaA_pkaC"]["all"]:
        tps = list(mean_distances["acaA_pkaC"]["all"].keys())
        ax.plot(
            tps,
            list(mean_distances["acaA_pkaC"]["all"].values()),
            label="acaA− pkaCOE",
            color="black",
            linestyle="dotted",
            marker="o",
            markerfacecolor="none",
        )
        add_error_bars(ax, "acaA_pkaC", tps, "black")
    if mean_distances["acaA_pkaC_psp"]["all"]:
        tps = list(mean_distances["acaA_pkaC_psp"]["all"].keys())
        ax.plot(
            tps,
            list(mean_distances["acaA_pkaC_psp"]["all"].values()),
            label="acaA− pkaCOE prespore",
            color="deepskyblue",
            linestyle="dotted",
            marker="o",
            markerfacecolor="none",
        )
        add_error_bars(ax, "acaA_pkaC_psp", tps, "deepskyblue")
    if mean_distances["acaA_pkaC_pst"]["all"]:
        tps = list(mean_distances["acaA_pkaC_pst"]["all"].keys())
        ax.plot(
            tps,
            list(mean_distances["acaA_pkaC_pst"]["all"].values()),
            label="acaA− pkaCOE prestalk",
            color="red",
            linestyle="dotted",
            marker="o",
            markerfacecolor="none",
        )
        add_error_bars(ax, "acaA_pkaC_pst", tps, "red")
    ax.set_xlabel("Time (h)", fontsize=18)
    ax.set_ylabel("Synchronicity (arbitrary units)", fontsize=18)
    if CI_dict is None:
        ax.legend(fontsize=10)
    else:
        ax.legend(fontsize=14)
    # Format x-tick labels
    if mean_distances["AX4"]["all"]:
        tps = list(mean_distances["AX4"]["all"].keys())
        ax.set_xticks(range(len(tps)))
        ax.set_xticklabels([format_timepoint(tp) for tp in tps], fontsize=14)
    ax.tick_params(axis='y', labelsize=14)

    # Panel 2: AX4
    ax = axes[1]
    if mean_distances["AX4"]["all"]:
        tps = list(mean_distances["AX4"]["all"].keys())
        ax.plot(
            tps,
            list(mean_distances["AX4"]["all"].values()),
            label="all cells" if not CI_dict else "AX4",
            color="black",
            marker="o",
        )
        add_error_bars(ax, "AX4", tps, "black")
    for marker in ax4_markers:
        if marker in mean_distances["AX4"] and mean_distances["AX4"][marker]:
            if CI_dict is None:
                ax.plot(
                    list(mean_distances["AX4"][marker].keys()),
                    list(mean_distances["AX4"][marker].values()),
                    label=f"all {marker_labels.get(marker, marker)}",
                    color=marker_colors.get(marker, "gray"),
                    alpha=0.4,
                )
            else:
                ax.plot(
                    list(mean_distances["AX4"][marker].keys()),
                    list(mean_distances["AX4"][marker].values()),
                    color="black",
                    marker="o",
                    linestyle="None",
                    markeredgewidth=0,
                    alpha=0.4,
                )
    if mean_distances["AX4_psp"]["all"]:
        tps = list(mean_distances["AX4_psp"]["all"].keys())
        ax.plot(
            tps,
            list(mean_distances["AX4_psp"]["all"].values()),
            label="prespore cells" if not CI_dict else "AX4 prespore",
            color="deepskyblue",
            marker="o",
        )
        add_error_bars(ax, "AX4_psp", tps, "deepskyblue")
        for marker in ax4_markers:
            if marker in mean_distances["AX4_psp"]:
                if (
                    marker in mean_distances["AX4_psp"]
                    and mean_distances["AX4_psp"][marker]
                ):
                    if CI_dict is None:
                        ax.plot(
                            list(mean_distances["AX4_psp"][marker].keys()),
                            list(mean_distances["AX4_psp"][marker].values()),
                            label=f"prespore {marker_labels.get(marker, marker)}",
                            color=marker_colors.get(marker, "gray"),
                            alpha=0.4,
                            linestyle="dashed",
                        )
                    else:
                        ax.plot(
                            list(mean_distances["AX4_psp"][marker].keys()),
                            list(mean_distances["AX4_psp"][marker].values()),
                            color="deepskyblue",
                            marker="o",
                            linestyle="None",
                            markeredgewidth=0,
                            alpha=0.4,
                        )
    if mean_distances["AX4_pst"]["all"]:
        tps = list(mean_distances["AX4_pst"]["all"].keys())
        ax.plot(
            tps,
            list(mean_distances["AX4_pst"]["all"].values()),
            label="prestalk cells" if not CI_dict else "AX4 prestalk",
            color="red",
            marker="o",
        )
        add_error_bars(ax, "AX4_pst", tps, "red")
        for marker in ax4_markers:
            if marker in mean_distances["AX4_pst"]:
                if (
                    marker in mean_distances["AX4_pst"]
                    and mean_distances["AX4_pst"][marker]
                ):
                    if CI_dict is None:
                        ax.plot(
                            list(mean_distances["AX4_pst"][marker].keys()),
                            list(mean_distances["AX4_pst"][marker].values()),
                            label=f"prestalk {marker_labels.get(marker, marker)}",
                            color=marker_colors.get(marker, "gray"),
                            alpha=0.4,
                            linestyle="dotted",
                        )
                    else:
                        ax.plot(
                            list(mean_distances["AX4_pst"][marker].keys()),
                            list(mean_distances["AX4_pst"][marker].values()),
                            color="red",
                            marker="o",
                            linestyle="None",
                            markeredgewidth=0,
                            alpha=0.4,
                        )
    ax.set_xlabel("Time (h)", fontsize=18)
    ax.set_ylabel("Synchronicity (arbitrary units)", fontsize=18)
    ax.set_title("AX4", fontsize=21)
    if CI_dict is None:
        ax.legend(fontsize=10)
    else:
        ax.legend(fontsize=14)
    # Format x-tick labels
    if mean_distances["AX4"]["all"]:
        tps = list(mean_distances["AX4"]["all"].keys())
        ax.set_xticks(range(len(tps)))
        ax.set_xticklabels([format_timepoint(tp) for tp in tps], fontsize=14)
    ax.tick_params(axis='y', labelsize=14)

    # Panel 3: acaA
    ax = axes[2]
    if mean_distances["acaA"]["all"]:
        tps = list(mean_distances["acaA"]["all"].keys())
        ax.plot(
            tps,
            list(mean_distances["acaA"]["all"].values()),
            label="all cells" if not CI_dict else "acaA−",
            color="black",
            marker="o",
            linestyle="solid" if not CI_dict else "dashed",
        )
        add_error_bars(ax, "acaA", tps, "black")
    for marker in acaA_markers:
        if marker in mean_distances["acaA"] and mean_distances["acaA"][marker]:
            if CI_dict is None:
                ax.plot(
                    list(mean_distances["acaA"][marker].keys()),
                    list(mean_distances["acaA"][marker].values()),
                    label=f"all {marker_labels.get(marker, marker)}",
                    color=marker_colors.get(marker, "gray"),
                    alpha=0.4,
                )
            else:
                ax.plot(
                    list(mean_distances["acaA"][marker].keys()),
                    list(mean_distances["acaA"][marker].values()),
                    color="black",
                    marker="o",
                    linestyle="None",
                    markeredgewidth=0,
                    alpha=0.4,
                )
    ax.set_xlabel("Time (h)", fontsize=18)
    ax.set_ylabel("Synchronicity (arbitrary units)", fontsize=18)
    ax.set_title("acaA−", fontsize=21)
    if CI_dict is None:
        ax.legend(fontsize=10)
    else:
        ax.legend(fontsize=14)
    # Format x-tick labels
    if mean_distances["acaA"]["all"]:
        tps = list(mean_distances["acaA"]["all"].keys())
        ax.set_xticks(range(len(tps)))
        ax.set_xticklabels([format_timepoint(tp) for tp in tps], fontsize=14)
    ax.tick_params(axis='y', labelsize=14)

    # Panel 4: acaA_pkaC
    ax = axes[3]
    if mean_distances["acaA_pkaC"]["all"]:
        tps = list(mean_distances["acaA_pkaC"]["all"].keys())
        ax.plot(
            tps,
            list(mean_distances["acaA_pkaC"]["all"].values()),
            label="all cells" if not CI_dict else "acaA− pkaCOE",
            color="black",
            marker="o",
            markerfacecolor="none" if CI_dict else "black",
            linestyle="dotted" if CI_dict else "solid",
        )
        add_error_bars(ax, "acaA_pkaC", tps, "black")
    for marker in acaA_pkaC_markers:
        if (
            marker in mean_distances["acaA_pkaC"]
            and mean_distances["acaA_pkaC"][marker]
        ):
            if CI_dict is None:
                ax.plot(
                    list(mean_distances["acaA_pkaC"][marker].keys()),
                    list(mean_distances["acaA_pkaC"][marker].values()),
                    label=f"all {marker_labels.get(marker, marker)}",
                    color=marker_colors.get(marker, "gray"),
                    alpha=0.4,
                )
            else:
                ax.plot(
                    list(mean_distances["acaA_pkaC"][marker].keys()),
                    list(mean_distances["acaA_pkaC"][marker].values()),
                    color="black",
                    marker="o",
                    linestyle="None",
                    markerfacecolor="none",
                    alpha=0.4,
                )

    if mean_distances["acaA_pkaC_psp"]["all"]:
        tps = list(mean_distances["acaA_pkaC_psp"]["all"].keys())
        ax.plot(
            tps,
            list(mean_distances["acaA_pkaC_psp"]["all"].values()),
            label="prespore cells" if not CI_dict else "acaA− pkaCOE prespore",
            color="deepskyblue",
            marker="o",
            markerfacecolor="none" if CI_dict else "deepskyblue",
            linestyle="dotted" if CI_dict else "solid",
        )
        add_error_bars(ax, "acaA_pkaC_psp", tps, "deepskyblue")
        for marker in acaA_pkaC_markers:
            if (
                marker in mean_distances["acaA_pkaC_psp"]
                and mean_distances["acaA_pkaC_psp"][marker]
            ):
                if CI_dict is None:
                    ax.plot(
                        list(mean_distances["acaA_pkaC_psp"][marker].keys()),
                        list(mean_distances["acaA_pkaC_psp"][marker].values()),
                        label=f"prespore {marker_labels.get(marker, marker)}",
                        color=marker_colors.get(marker, "gray"),
                        alpha=0.4,
                        linestyle="dashed",
                    )
                else:
                    ax.plot(
                        list(mean_distances["acaA_pkaC_psp"][marker].keys()),
                        list(mean_distances["acaA_pkaC_psp"][marker].values()),
                        color="deepskyblue",
                        marker="o",
                        linestyle="None",
                        markerfacecolor="none",
                        alpha=0.4,
                    )
    if mean_distances["acaA_pkaC_pst"]["all"]:
        tps = list(mean_distances["acaA_pkaC_pst"]["all"].keys())
        ax.plot(
            tps,
            list(mean_distances["acaA_pkaC_pst"]["all"].values()),
            label="prestalk cells" if not CI_dict else "acaA− pkaCOE prestalk",
            color="red",
            marker="o",
            markerfacecolor="none" if CI_dict else "red",
            linestyle="dotted" if CI_dict else "solid",
        )
        add_error_bars(ax, "acaA_pkaC_pst", tps, "red")
        for marker in acaA_pkaC_markers:
            if (
                marker in mean_distances["acaA_pkaC_pst"]
                and mean_distances["acaA_pkaC_pst"][marker]
            ):
                if CI_dict is None:
                    ax.plot(
                        list(mean_distances["acaA_pkaC_pst"][marker].keys()),
                        list(mean_distances["acaA_pkaC_pst"][marker].values()),
                        label=f"prestalk {marker_labels.get(marker, marker)}",
                        color=marker_colors.get(marker, "gray"),
                        alpha=0.4,
                        linestyle="dotted",
                    )
                else:
                    ax.plot(
                        list(mean_distances["acaA_pkaC_pst"][marker].keys()),
                        list(mean_distances["acaA_pkaC_pst"][marker].values()),
                        color="red",
                        marker="o",
                        linestyle="None",
                        markerfacecolor="none",
                        alpha=0.4,
                    )
    ax.set_xlabel("Time (h)", fontsize=18)
    ax.set_ylabel("Synchronicity (arbitrary units)", fontsize=18)
    ax.set_title("acaA− pkaCOE", fontsize=21)
    if CI_dict is None:
        ax.legend(fontsize=10)
    else:
        ax.legend(fontsize=14)
    # Format x-tick labels
    if mean_distances["acaA_pkaC"]["all"]:
        tps = list(mean_distances["acaA_pkaC"]["all"].keys())
        ax.set_xticks(range(len(tps)))
        ax.set_xticklabels([format_timepoint(tp) for tp in tps], fontsize=14)
    ax.tick_params(axis='y', labelsize=14)

    if suptitle:
        fig.suptitle(suptitle, fontsize=14)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        plt.savefig(save_path.replace(".png", ".pdf"), dpi=300, bbox_inches="tight")

    plt.show()
    return fig, axes


def export_synchronicity_data(
    results: Dict,
    CI_dict: Optional[Dict] = None,
    save_path: str = "synchronicity_source_data.xlsx",
):
    """
    Export synchronicity data to Excel format for publication.
    
    Parameters
    ----------
    results : dict
        Output from compute_synchronicity() containing:
        - 'mean_distances': synchronicity values
        - 'cell_counts': cell counts
        - 'median_genes': median genes per cell
        - 'total_reads': total reads
    CI_dict : dict, optional
        Output from compute_confidence_intervals()
    save_path : str
        Path to save the Excel file
    
    Returns
    -------
    dict of DataFrames
        Dictionary containing all the DataFrames that were saved
    """
    
    mean_distances = results['mean_distances']
    
    dfs = {}
    
    # Helper to convert nested dict to DataFrame
    def dict_to_df(data_dict, value_name="value"):
        rows = []
        for strain, subsets in data_dict.items():
            for subset, timepoints in subsets.items():
                for tp, value in timepoints.items():
                    rows.append({
                        "strain": strain,
                        "subset": subset,
                        "timepoint": tp,
                        value_name: value
                    })
        return pd.DataFrame(rows)
    
    # Create DataFrames for each metric
    dfs['synchronicity'] = dict_to_df(mean_distances, "synchronicity")
    
    # Add confidence intervals if provided
    if CI_dict is not None:
        ci_rows = []
        for strain, timepoints in CI_dict.items():
            for tp, values in timepoints.items():
                if len(values) >= 3:
                    # Get the actual plotted synchronicity value
                    sync_value = mean_distances.get(strain, {}).get("all", {}).get(tp, None)
                    if sync_value is not None:
                        # Error bar distances (same as in plot)
                        error_lower = values[0] - values[1]  # mean - lower
                        error_upper = values[2] - values[0]  # upper - mean
                        # Actual plotted CI bounds
                        ci_lower_plotted = sync_value - error_lower
                        ci_upper_plotted = sync_value + error_upper
                        ci_rows.append({
                            "strain": strain,
                            "timepoint": tp,
                            "synchronicity": sync_value,
                            "CI_lower": ci_lower_plotted,
                            "CI_upper": ci_upper_plotted,
                        })
        dfs['synchronicity_confidence_intervals'] = pd.DataFrame(ci_rows)
    
    # Save to Excel with multiple sheets
    with pd.ExcelWriter(save_path, engine='openpyxl') as writer:
        for sheet_name, df in dfs.items():
            # Excel sheet names have max 31 characters
            sheet_name_truncated = sheet_name[:31]
            df.to_excel(writer, sheet_name=sheet_name_truncated, index=True)
    
    print(f"Data exported to {save_path}")
    print(f"Sheets: {list(dfs.keys())}")
    
    return dfs
