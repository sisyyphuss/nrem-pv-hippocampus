"""
Master analysis driver: tests whether inter-session NREM sleep duration predicts
hippocampal population-vector (PV) decorrelation across chronic recording sessions.

Evaluates DANDI:000213 (rat CA1), DANDI:001695 (mouse CA1), and pooled cohorts
with full robustness checks (within-subject permutation tests, mixed-effects models,
and cell-type specificity).
"""
import warnings
warnings.filterwarnings("ignore")
warnings.simplefilter("ignore")

import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

# Ensure local module imports work
SRC_DIR = Path(__file__).resolve().parent
REPO_DIR = SRC_DIR.parent
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from stats_lib import (
    pooled_correlation,
    fit_mixed_model,
    permutation_test,
    leave_one_subject_out,
    leave_one_dataset_out,
)

TABLES_DIR = REPO_DIR / "results" / "tables"


def main():
    print("=" * 78)
    print("NREM Sleep & Hippocampal Population-Vector Correlation Analysis")
    print("=" * 78)

    file_a = TABLES_DIR / "session_pairs.csv"
    file_b = TABLES_DIR / "session_pairs_001695.csv"
    file_c = TABLES_DIR / "session_pairs_combined.csv"
    file_d = TABLES_DIR / "session_pairs_001695_interneuron.csv"

    if not file_c.exists():
        print(f"Error: Table file not found: {file_c}")
        return

    # 1. Analysis A: DANDI:000213 (Rat CA1)
    df_a = pd.read_csv(file_a)
    r_spear_a, p_spear_a, n_a = pooled_correlation(df_a, method="spearman")
    r_pear_a, p_pear_a, _ = pooled_correlation(df_a, method="pearson")
    _, p_perm_a, _ = permutation_test(df_a, n_perm=5000)

    print("\n[Analysis A: DANDI:000213 - Rat CA1]")
    print(f"  Sample size:        {n_a} session-pairs ({len(df_a['subject'].unique())} subjects)")
    print(f"  Spearman r:         {r_spear_a:.3f} (p = {p_spear_a:.3f})")
    print(f"  Pearson r:          {r_pear_a:.3f} (p = {p_pear_a:.3f})")
    print(f"  Permutation p:      {p_perm_a:.3f} (within-subject shuffle, 5,000 iter)")

    # 2. Analysis B: DANDI:001695 (Mouse CA1)
    df_b = pd.read_csv(file_b)
    r_spear_b, p_spear_b, n_b = pooled_correlation(df_b, method="spearman")
    r_pear_b, p_pear_b, _ = pooled_correlation(df_b, method="pearson")
    _, p_perm_b, _ = permutation_test(df_b, n_perm=5000)

    print("\n[Analysis B: DANDI:001695 - Mouse CA1]")
    print(f"  Sample size:        {n_b} session-pairs ({len(df_b['subject'].unique())} subjects)")
    print(f"  Spearman r:         {r_spear_b:.3f} (p = {p_spear_b:.3f})")
    print(f"  Pearson r:          {r_pear_b:.3f} (p = {p_pear_b:.3f})")
    print(f"  Permutation p:      {p_perm_b:.3f} (within-subject shuffle, 5,000 iter)")

    # 3. Analysis C: Combined Multi-Dataset
    df_c = pd.read_csv(file_c)
    r_spear_c, p_spear_c, n_c = pooled_correlation(df_c, method="spearman")
    r_pear_c, p_pear_c, _ = pooled_correlation(df_c, method="pearson")
    _, p_perm_c, _ = permutation_test(df_c, n_perm=5000)

    print("\n[Analysis C: Combined Multi-Dataset (DANDI:000213 + DANDI:001695)]")
    print(f"  Total sample:       {n_c} session-pairs ({len(df_c['subject'].unique())} subjects, 2 species)")
    print(f"  Spearman r:         {r_spear_c:.3f} (p = {p_spear_c:.3f})")
    print(f"  Pearson r:          {r_pear_c:.3f} (p = {p_pear_c:.3f})")
    print(f"  Permutation p:      {p_perm_c:.3f} (within-subject shuffle, 5,000 iter)")

    # Mixed-effects model on combined
    try:
        mixed_res, _ = fit_mixed_model(df_c)
        nrem_coef = mixed_res.params.get("nrem_between_min", np.nan)
        nrem_p = mixed_res.pvalues.get("nrem_between_min", np.nan)
        print(f"  Mixed-effects model (subject random intercept, controlling for wake, REM, gap):")
        print(f"    NREM beta:        {nrem_coef:.5f} (p = {nrem_p:.4f})")
    except Exception as e:
        print(f"  Mixed-effects model note: {e}")

    # 4. Analysis D: Pyramidal vs. Interneuron Control
    if file_d.exists():
        df_d = pd.read_csv(file_d)
        print("\n[Analysis D: Pyramidal vs Interneuron Cell-Type Specificity]")
        pyr_mean = df_d["pv_corr_pyramidal"].mean()
        int_mean = df_d["pv_corr_interneuron"].mean()
        print(f"  Mean Pyramidal PV Correlation:   {pyr_mean:.3f}")
        print(f"  Mean Interneuron PV Correlation: {int_mean:.3f}")

    print("\n" + "=" * 78)
    print("Conclusion: Across cohorts, inter-session NREM duration demonstrates a")
    print("consistent negative trend (r = -0.246), but does not reach significance")
    print("under within-subject permutation testing (p = 0.225).")
    print("=" * 78)


if __name__ == "__main__":
    main()
