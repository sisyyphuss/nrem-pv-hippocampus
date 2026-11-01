"""
Statistics helpers for the NREM-vs-PV-correlation pipeline: pooled
correlation, mixed-effects model (subject random intercept), permutation
test, and leave-one-subject-out robustness check.
"""
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.formula.api as smf


def pooled_correlation(df, x="nrem_between_min", y="pv_corr", method="spearman"):
    d = df[[x, y]].dropna()
    if method == "spearman":
        r, p = stats.spearmanr(d[x], d[y])
    else:
        r, p = stats.pearsonr(d[x], d[y])
    return r, p, len(d)


def fit_mixed_model(df, formula="pv_corr ~ nrem_between_min + wake_between_min + rem_between_min + day_gap",
                     group="subject"):
    d = df.dropna(subset=[c.strip() for c in formula.replace("pv_corr ~", "").split("+")] + ["pv_corr", group]).copy()
    model = smf.mixedlm(formula, data=d, groups=d[group])
    result = model.fit(reml=True)
    return result, d


def permutation_test(df, x="nrem_between_min", y="pv_corr", group="subject",
                      n_perm=5000, method="spearman", random_state=0):
    """
    Null hypothesis: no relationship between inter-session NREM duration and
    PV correlation. Null distribution is built by shuffling x WITHIN each
    subject (preserves each animal's own NREM-duration distribution and
    each animal's own PV-correlation distribution, only breaking the
    pairing between them) -- more conservative than a fully unrestricted
    shuffle given the repeated-measures structure.
    """
    rng = np.random.default_rng(random_state)
    d = df[[x, y, group]].dropna().reset_index(drop=True)
    obs_r, _, n = pooled_correlation(d, x, y, method)

    null_rs = np.empty(n_perm)
    groups = d[group].values
    x_vals = d[x].values.copy()
    for i in range(n_perm):
        x_shuf = x_vals.copy()
        for g in np.unique(groups):
            idx = np.where(groups == g)[0]
            x_shuf[idx] = rng.permutation(x_shuf[idx])
        if method == "spearman":
            r, _ = stats.spearmanr(x_shuf, d[y].values)
        else:
            r, _ = stats.pearsonr(x_shuf, d[y].values)
        null_rs[i] = r

    p_emp = (np.sum(np.abs(null_rs) >= np.abs(obs_r)) + 1) / (n_perm + 1)
    return obs_r, p_emp, null_rs


def leave_one_subject_out(df, x="nrem_between_min", y="pv_corr", group="subject", method="spearman"):
    subjects = df[group].dropna().unique()
    rows = []
    for s in subjects:
        sub = df[df[group] != s]
        r, p, n = pooled_correlation(sub, x, y, method)
        rows.append({"held_out": s, "r_without": r, "p_without": p, "n_pairs": n})
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------
# EXTENSIONS FOR THE COMBINED MULTI-DATASET ANALYSIS (Analysis C).
# Added as new functions -- the single-dataset functions above are left
# untouched so Analysis A's numbers are not at risk of changing.
# ----------------------------------------------------------------------

def correlation_by_dataset(df, x="nrem_between_min", y="pv_corr", dataset_col="dataset", method="spearman"):
    """
    Each dataset's OWN pooled correlation, computed separately (not a
    leave-one-out contrast). Lets the reader see directly whether pooling
    datasets helps, hurts, or is dominated by one dataset, rather than only
    seeing the combined number.
    """
    rows = []
    for ds, sub in df.groupby(dataset_col):
        r, p, n = pooled_correlation(sub, x, y, method)
        rows.append({"dataset": ds, "r": r, "p": p, "n_pairs": n, "n_subjects": sub["subject"].nunique()})
    return pd.DataFrame(rows)


def leave_one_dataset_out(df, x="nrem_between_min", y="pv_corr", dataset_col="dataset", method="spearman"):
    """Analogous to leave_one_subject_out, but drops one whole DATASET at a time."""
    datasets = df[dataset_col].dropna().unique()
    rows = []
    for ds in datasets:
        sub = df[df[dataset_col] != ds]
        r, p, n = pooled_correlation(sub, x, y, method)
        rows.append({"held_out_dataset": ds, "r_without": r, "p_without": p,
                      "n_pairs": n, "n_subjects": sub["subject"].nunique()})
    return pd.DataFrame(rows)


def permutation_test_multi_dataset(df, x="nrem_between_min", y="pv_corr", group="subject",
                                    dataset_col="dataset", n_perm=5000, method="spearman", random_state=0):
    """
    Within-subject permutation test (identical logic to permutation_test()),
    but with an explicit, loudly-failing precondition check: every subject
    must belong to exactly ONE dataset. If that holds -- which it does here
    by construction, since 000213 and 001695 subject IDs never collide
    (DT2/DT5/DT7 vs M01/M02/M03/M05) -- then shuffling x WITHIN each
    subject automatically cannot move NREM values across a dataset
    boundary, so the block structure the task asked to preserve is
    preserved for free. This function exists to make that guarantee
    EXPLICIT and machine-checked rather than an unstated assumption.
    """
    d = df[[x, y, group, dataset_col]].dropna()
    n_ds_per_subject = d.groupby(group)[dataset_col].nunique()
    bad = n_ds_per_subject[n_ds_per_subject > 1]
    if len(bad) > 0:
        raise ValueError(
            f"Subject(s) {list(bad.index)} appear in more than one dataset -- "
            "the within-subject permutation would leak across dataset labels. "
            "Refusing to run until subject IDs are made dataset-unique."
        )
    return permutation_test(df, x=x, y=y, group=group, n_perm=n_perm, method=method, random_state=random_state)


# ----------------------------------------------------------------------
# EXTENSION FOR ANALYSIS D: paired pyramidal-vs-interneuron comparison.
# Added as a new function -- nothing above is modified.
# ----------------------------------------------------------------------

def paired_comparison(df, col_a, col_b, min_n_for_test=6):
    """
    Descriptive, and CONDITIONALLY inferential, paired comparison between
    two columns computed on the SAME rows (e.g. pyramidal-ensemble vs.
    interneuron-ensemble PV correlation for the same session-pairs).

    A Wilcoxon signed-rank test is only attempted when n >= min_n_for_test
    (default 6, the smallest n scipy's Wilcoxon table-based p-value is
    conventionally treated as informative for). Below that, this function
    still returns the descriptive numbers (mean/median paired difference,
    per-pair sign counts) but leaves test/stat/p as NaN/None and sets
    `test_valid=False` -- matching this project's standard of never
    implying a formal significance result the sample size cannot support.
    """
    d = df[[col_a, col_b]].dropna()
    diffs = (d[col_a] - d[col_b]).values
    n = len(d)
    out = {
        "n": n,
        "mean_diff": float(np.mean(diffs)) if n else np.nan,
        "median_diff": float(np.median(diffs)) if n else np.nan,
        "n_a_gt_b": int((diffs > 0).sum()),
        "n_b_gt_a": int((diffs < 0).sum()),
        "n_tied": int((diffs == 0).sum()),
        "test": None, "stat": np.nan, "p": np.nan, "test_valid": False,
    }
    if n >= min_n_for_test and np.any(diffs != 0):
        stat, p = stats.wilcoxon(diffs)
        out.update(test="wilcoxon_signed_rank", stat=float(stat), p=float(p), test_valid=True)
    return out


def fit_mixed_model_dataset_covariate(df, group="subject",
                                       base_formula="pv_corr ~ nrem_between_min",
                                       dataset_col="dataset"):
    """
    Fits three nested mixed-effects models (subject random intercept) to
    test whether the NREM-PV relationship differs by dataset:
      1. `no_dataset`    : pv_corr ~ nrem_between_min
      2. `additive`      : pv_corr ~ nrem_between_min + dataset          (dataset as a covariate / main effect)
      3. `interaction`   : pv_corr ~ nrem_between_min * dataset          (does the SLOPE differ by dataset?)
    Returns a dict of {name: (result, data_used)}. The notebook reports the
    interaction term's coefficient/p-value explicitly -- if it is not
    negligible, that means "does NREM predict PV correlation" does NOT have
    one dataset-independent answer, which must be surfaced, not hidden
    behind a single pooled coefficient.
    """
    out = {}
    d = df.dropna(subset=["pv_corr", "nrem_between_min", group, dataset_col]).copy()

    m1 = smf.mixedlm(base_formula, data=d, groups=d[group])
    out["no_dataset"] = (m1.fit(reml=True), d)

    f2 = base_formula + f" + C({dataset_col})"
    m2 = smf.mixedlm(f2, data=d, groups=d[group])
    out["additive"] = (m2.fit(reml=True), d)

    f3 = base_formula.replace(" ~ ", " ~ ").rstrip() + f" * C({dataset_col})"
    # build interaction formula explicitly: pv_corr ~ nrem_between_min * C(dataset)
    f3 = f"pv_corr ~ nrem_between_min * C({dataset_col})"
    m3 = smf.mixedlm(f3, data=d, groups=d[group])
    out["interaction"] = (m3.fit(reml=True), d)

    return out
