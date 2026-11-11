"""Histogram helpers shared by every analysis script (and the tests).

These used to live in analysis/m01_analysis_uproot.py and were imported from
there by the closure test; the project rule is that anything used twice lives
in src/hicpy/.

All helpers cast their inputs to float64 first: ROOT branches are float32, and
numpy's weighted histogram accumulates a cumulative sum in the weights' dtype,
which silently quantises sparsely populated bins (docs/troubleshooting.md #12).
"""
from __future__ import annotations

import numpy as np


def make_edges(spec):
    """Bin edges from a config block {type: log|linear, n_bins, lo, hi}."""
    if spec["type"] == "log":
        return np.geomspace(spec["lo"], spec["hi"], spec["n_bins"] + 1)
    return np.linspace(spec["lo"], spec["hi"], spec["n_bins"] + 1)


def centers_widths(edges):
    """Arithmetic bin centres and bin widths."""
    edges = np.asarray(edges, dtype=np.float64)
    return 0.5 * (edges[1:] + edges[:-1]), np.diff(edges)


def hist_counts(values, edges):
    """Unweighted histogram -> (counts, sqrt(N) errors). Poisson errors are the
    correct statistical uncertainty for a raw count in a bin."""
    counts, _ = np.histogram(np.asarray(values, dtype=np.float64), bins=edges)
    counts = counts.astype(np.float64)
    return counts, np.sqrt(counts)


def hist_weighted(values, weights, edges):
    """Weighted histogram -> (sum of weights, sqrt(sum of weights^2)) per bin.

    This is exactly what ROOT's TH1::Sumw2() does.  With per-particle weights
    the Poisson sqrt(N) rule no longer applies: the correct statistical error on
    a weighted sum is sqrt(sum w_i^2).
    """
    values = np.asarray(values, dtype=np.float64)      # troubleshooting #12
    weights = np.asarray(weights, dtype=np.float64)
    sumw, _ = np.histogram(values, bins=edges, weights=weights)
    sumw2, _ = np.histogram(values, bins=edges, weights=weights ** 2)
    return sumw, np.sqrt(sumw2)


def hist_record(edges, contents, errors):
    """Serialisable histogram (edges, per-bin contents and errors).

    Every arm writes its histograms in this form so that m01_compare_arms.py
    can compare them bin by bin -- agreeing on a handful of scalar summaries is
    not evidence that two implementations agree on the histograms.
    """
    return {"edges": [float(v) for v in edges],
            "contents": [float(v) for v in contents],
            "errors": [float(v) for v in errors]}


def select_bins(edges, lo, hi):
    """Boolean mask of the bins whose CENTRE lies strictly inside (lo, hi).

    The same rule as ROOT's TH1::Fit "R" option, so that every arm and every
    fit method uses exactly the same bins for a stated fit range.
    """
    ctr, _ = centers_widths(edges)
    return (ctr > lo) & (ctr < hi)


# Backwards-compatible name used by the Session 1-3 scripts.
hist_with_errors = hist_counts

__all__ = ["make_edges", "centers_widths", "hist_counts", "hist_with_errors",
           "hist_weighted", "hist_record", "select_bins"]
