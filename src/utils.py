import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import gaussian_kde
from scipy.signal import argrelextrema
from sklearn.metrics.pairwise import euclidean_distances


def assign_markers(adata, time_points, neon_all=False, verbose=False):
    """
    For each hr in time_points, computes per-marker density-based thresholds,
    aligns them, assigns each cell its top marker or 'Unassigned', and stores
    the result in adata.obs['marker'].

    Parameters
    ----------
    adatas : dict of AnnData
    time_points : iterable of str
    neon_all : bool
        If True, always include mNeonG; otherwise skip it for '00hr' and '16hr'.
    verbose : bool
        If True, prints hours, min cell counts, thresholds, and final counts.
    """
    adatas = {}
    for hr in time_points:
        adatas[hr] = adata[adata.obs["time"] == hr]
        obs = adatas[hr].obs
        # pick genes
        genes = ["act15GFP", "mCherry", "mCerulean"]
        if neon_all or hr not in ("00hr", "16hr"):
            genes.append("mNeonG")
        # build expression frame
        df = pd.DataFrame({g: np.log1p(obs[f"{g}"]) for g in genes})
        # count non‐zero cells & report
        nz = {g: (df[g] > 0).sum() for g in genes}
        if verbose:
            print(f"Hours: {hr} — min cells: {min(nz.values())}")
        # find thresholds
        th = {}
        for g in genes:
            x = df[g]
            x = x[x > 0]
            if x.empty:
                th[g] = 0
            else:
                kde = gaussian_kde(x)
                grid = np.linspace(x.min(), x.max(), 1_000)
                y = kde(grid)
                peaks = argrelextrema(y, np.greater)[0]
                if peaks.size:
                    # pick first non‐background peak
                    t0 = grid[peaks[0]]
                    if t0 < 0.2 and peaks.size > 1:
                        t0 = grid[peaks[1]]
                    th[g] = t0
                else:
                    th[g] = np.percentile(x, 99)
            if verbose:
                print(f"  {g} threshold = {th[g]:.3f}")
        # align thresholds
        m0 = min(th.values())
        adj = {g: th[g] - m0 for g in genes}
        # score & assign
        scores = pd.DataFrame({g: df[g].sub(adj[g]) for g in genes}, index=obs.index)
        best = scores.idxmax(axis=1)
        mask = scores.max(axis=1) > 0
        adatas[hr].obs["marker"] = np.where(mask, best, "Unassigned")
        if verbose:
            print(adatas[hr].obs["marker"].value_counts(), "\n")
    # combine
    adata_all = sc.concat(adatas, join="outer")
    adata.obs["marker"] = adata_all.obs["marker"]

    return adata


def average_pairwise_euclidean_similarity(X):
    """
    Computes the average pairwise euclidean similarity among all rows of X.

    Parameters
    ----------
    X : np.ndarray
        N x D matrix where each row is a cell embedding.

    Returns
    -------
    float
        The average pairwise euclidean similarity.
    """
    n = X.shape[0]
    if n < 2:
        return np.nan
    S = euclidean_distances(X)
    S_off_diag = S[np.triu_indices(n, k=1)]
    avg_similarity = (2 - S_off_diag.mean()) / 2
    return avg_similarity

def bootstrap_confidence_interval(data, num_bootstrap=10000, ci=95):
    """
    Calculate bootstrap confidence intervals for the mean of the data.

    Args:
        data: Array-like, the data to bootstrap.
        num_bootstrap: Number of bootstrap samples.
        ci: Confidence interval percentage (default is 95%).

    Returns:
        mean: Mean of the original data.
        lower: Lower bound of the confidence interval.
        upper: Upper bound of the confidence interval.
    """
    bootstrap_means = []
    n = len(data)
    for _ in range(num_bootstrap):
        sample = np.random.choice(data, size=n, replace=True)  # Resample with replacement
        bootstrap_means.append(np.mean(sample))
    
    lower = np.percentile(bootstrap_means, (100 - ci) / 2)
    upper = np.percentile(bootstrap_means, 100 - (100 - ci) / 2)
    return [np.mean(data), lower, upper, bootstrap_means]