"""
Figure-generation functions for the NREM-vs-PV-correlation notebook.
Matplotlib only, no seaborn dependency, colorblind-safe palette,
dark/light-mode-agnostic (explicit colors, not rcParam defaults).
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

# Palette: colorblind-safe categorical set (Okabe-Ito), consistent across figures
COLORS = {
    "DT2": "#0072B2",
    "DT5": "#D55E00",
    "DT7": "#009E73",
    "accent": "#CC79A7",
    "null": "#999999",
    "obs": "#D55E00",
    # DANDI:001695 subjects (mouse) -- kept visually distinct from the
    # DANDI:000213 rat subjects above by using the remaining Okabe-Ito hues.
    "M01": "#E69F00",
    "M02": "#56B4E9",
    "M03": "#F0E442",
    "M05": "#000000",
}

DATASET_COLORS = {"000213": "#0072B2", "001695": "#E69F00"}
DATASET_MARKERS = {"000213": "o", "001695": "^"}


def color_for_subject(subject):
    """Fallback-safe color lookup so new subjects never crash a figure call."""
    return COLORS.get(subject, "#888888")


def fig_hypnogram(sess, ax=None, title=""):
    if ax is None:
        fig, ax = plt.subplots(figsize=(9, 1.8))
    state_y = {"Awake": 2, "Non-REM": 1, "REM": 0}
    state_color = {"Awake": "#0072B2", "Non-REM": "#D55E00", "REM": "#009E73"}
    for lab, s0, s1 in zip(sess["ss_label"], sess["ss_start"], sess["ss_stop"]):
        lab = str(lab)
        if lab not in state_y:
            continue
        ax.barh(state_y[lab], s1 - s0, left=s0, height=0.8, color=state_color[lab])
    for name, col in [("before trials", "k"), ("during trials", "red"), ("after trials", "k")]:
        idx = np.where(sess["ep_labels"] == name)[0]
        if len(idx):
            ax.axvline(sess["ep_start"][idx[0]], color="0.3", lw=0.8, ls="--")
    ax.set_yticks([0, 1, 2]); ax.set_yticklabels(["REM", "NREM", "Wake"])
    ax.set_xlabel("Time in session (s)")
    ax.set_title(title, fontsize=10, loc="left")
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    return ax


def fig_rate_maps(rate_maps, occ_mask, n_show=6, title=""):
    idx = [i for i in range(rate_maps.shape[0]) if np.nanmean(rate_maps[i][occ_mask]) > 0.1]
    idx = idx[:n_show]
    n = len(idx)
    fig, axes = plt.subplots(1, max(n, 1), figsize=(2.2 * max(n, 1), 2.4))
    if n == 1:
        axes = [axes]
    for ax, i in zip(axes, idx):
        rm = np.where(occ_mask, rate_maps[i], np.nan)
        im = ax.imshow(rm.T, origin="lower", cmap="viridis", aspect="auto")
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f"unit {i}", fontsize=8)
    fig.suptitle(title, fontsize=10)
    fig.tight_layout()
    return fig


def fig_ensemble_tuning(vec_a, vec_b, label_a, label_b, r_val, title=""):
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.2))
    axes[0].plot(vec_a, color=COLORS["DT2"], lw=1.3, label=label_a)
    axes[0].plot(vec_b, color=COLORS["DT5"], lw=1.3, label=label_b)
    axes[0].set_xlabel("Spatial bin (flattened)"); axes[0].set_ylabel("Ensemble z-rate")
    axes[0].legend(frameon=False, fontsize=8)
    axes[0].set_title("Ensemble tuning curves", fontsize=10)

    ok = np.isfinite(vec_a) & np.isfinite(vec_b)
    axes[1].scatter(vec_a[ok], vec_b[ok], s=14, color=COLORS["accent"], alpha=0.8, edgecolor="none")
    lims = [min(vec_a[ok].min(), vec_b[ok].min()), max(vec_a[ok].max(), vec_b[ok].max())]
    axes[1].plot(lims, lims, color="0.6", lw=0.8, ls="--")
    axes[1].set_xlabel(label_a); axes[1].set_ylabel(label_b)
    axes[1].set_title(f"PV correlation r={r_val:.2f}", fontsize=10)
    for ax in axes:
        for spine in ["top", "right"]:
            ax.spines[spine].set_visible(False)
    fig.suptitle(title, fontsize=11)
    fig.tight_layout()
    return fig


def fig_main_scatter(df, x="nrem_between_min", y="pv_corr", r=None, p=None, method="Spearman"):
    fig, ax = plt.subplots(figsize=(5.5, 4.8))
    for subj, sub in df.groupby("subject"):
        ax.scatter(sub[x], sub[y], s=48, color=COLORS.get(subj, "gray"),
                   edgecolor="white", linewidth=0.6, label=f"{subj} (n={len(sub)})", zorder=3)
    d = df[[x, y]].dropna()
    if len(d) >= 3:
        slope, intercept, *_ = stats.linregress(d[x], d[y])
        xs = np.linspace(d[x].min(), d[x].max(), 50)
        ax.plot(xs, slope * xs + intercept, color="0.25", lw=1.6, zorder=2)
    ax.set_xlabel("NREM sleep recorded between sessions (min)")
    ax.set_ylabel("Population Vector correlation\n(ensemble tuning-curve method)")
    title = "NREM between sessions vs. PV correlation"
    if r is not None:
        title += f"\n{method} r={r:.2f}" + (f", p={p:.3f}" if p is not None else "")
    ax.set_title(title, fontsize=10)
    ax.legend(frameon=False, fontsize=8, loc="best")
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    return fig


def fig_permutation_null(null_rs, obs_r, p_emp):
    fig, ax = plt.subplots(figsize=(5.5, 3.8))
    ax.hist(null_rs, bins=40, color=COLORS["null"], alpha=0.8, label="null (within-subject shuffle)")
    ax.axvline(obs_r, color=COLORS["obs"], lw=2.2, label=f"observed r={obs_r:.2f}")
    ax.set_xlabel("Spearman r (NREM-between vs. PV correlation)")
    ax.set_ylabel("Count (permutations)")
    ax.set_title(f"Permutation test, empirical p={p_emp:.4f}", fontsize=10)
    ax.legend(frameon=False, fontsize=8)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    return fig


def fig_loso(loso_df, obs_r):
    fig, ax = plt.subplots(figsize=(5.5, 3.4))
    ax.axhline(obs_r, color="0.5", lw=1, ls="--", label="full-sample r")
    ax.bar(loso_df["held_out"], loso_df["r_without"], color=COLORS["accent"])
    ax.set_ylabel("r without this subject")
    ax.set_xlabel("Subject held out")
    ax.set_title("Leave-one-subject-out stability", fontsize=10)
    ax.legend(frameon=False, fontsize=8)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    return fig


def fig_method_robustness(df):
    fig, ax = plt.subplots(figsize=(5, 4.6))
    ok = df[["pv_corr", "pv_corr_rv"]].dropna()
    ax.scatter(ok["pv_corr"], ok["pv_corr_rv"], s=44, color=COLORS["DT7"], edgecolor="white", linewidth=0.5)
    if len(ok) >= 3:
        r, p = stats.spearmanr(ok["pv_corr"], ok["pv_corr_rv"])
        ax.set_title(f"Ensemble-tuning-curve r vs. RV-coefficient method\nSpearman r={r:.2f}, p={p:.3f}", fontsize=9)
    ax.set_xlabel("PV correlation (ensemble tuning-curve method)")
    ax.set_ylabel("RV coefficient (multivariate method)")
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    return fig


# ============================================================================
# EXTENSIONS FOR THE COMBINED MULTI-DATASET ANALYSIS (Analysis B / C) and the
# DANDI:000061 descriptive-only section. Added as new functions; nothing
# above this line was modified.
# ============================================================================

def fig_main_scatter_by_dataset(df, x="nrem_between_min", y="pv_corr", r=None, p=None, method="Spearman"):
    """Same idea as fig_main_scatter(), but colored/marker-coded by dataset
    (with subject-level color within each dataset) so it is immediately
    visible whether one dataset dominates the combined scatter."""
    fig, ax = plt.subplots(figsize=(6.2, 5.0))
    for (ds, subj), sub in df.groupby(["dataset", "subject"]):
        ax.scatter(sub[x], sub[y], s=50, color=color_for_subject(subj),
                   marker=DATASET_MARKERS.get(ds, "o"), edgecolor="white", linewidth=0.6,
                   label=f"{ds}:{subj} (n={len(sub)})", zorder=3)
    d = df[[x, y]].dropna()
    if len(d) >= 3:
        slope, intercept, *_ = stats.linregress(d[x], d[y])
        xs = np.linspace(d[x].min(), d[x].max(), 50)
        ax.plot(xs, slope * xs + intercept, color="0.25", lw=1.6, zorder=2)
    ax.set_xlabel("NREM sleep recorded between sessions (min)")
    ax.set_ylabel("Population Vector correlation\n(ensemble tuning-curve method)")
    title = "Combined (000213 + 001695): NREM between sessions vs. PV correlation"
    if r is not None:
        title += f"\n{method} r={r:.2f}" + (f", p={p:.3f}" if p is not None else "")
    ax.set_title(title, fontsize=10)
    ax.legend(frameon=False, fontsize=7, loc="best", ncol=2)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    return fig


def fig_correlation_by_dataset(corr_by_dataset_df, combined_r=None, combined_label="combined"):
    """Bar chart: each dataset's own pooled r, plus the combined r for comparison."""
    fig, ax = plt.subplots(figsize=(5.5, 3.8))
    xs = list(corr_by_dataset_df["dataset"])
    ys = list(corr_by_dataset_df["r"])
    colors = [DATASET_COLORS.get(d, "#888888") for d in xs]
    if combined_r is not None:
        xs = xs + [combined_label]
        ys = ys + [combined_r]
        colors = colors + ["0.3"]
    ax.bar(xs, ys, color=colors)
    ax.axhline(0, color="0.6", lw=0.8)
    ax.set_ylabel("Spearman r (NREM-between vs. PV correlation)")
    ax.set_title("Pooled correlation: each dataset alone vs. combined", fontsize=10)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    return fig


def fig_loso_dataset(lodo_df, obs_r):
    fig, ax = plt.subplots(figsize=(5, 3.4))
    ax.axhline(obs_r, color="0.5", lw=1, ls="--", label="full-sample r")
    ax.bar(lodo_df["held_out_dataset"], lodo_df["r_without"], color=COLORS["accent"])
    ax.set_ylabel("r without this dataset")
    ax.set_xlabel("Dataset held out")
    ax.set_title("Leave-one-dataset-out stability", fontsize=10)
    ax.legend(frameon=False, fontsize=8)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    return fig


def fig_pyr_vs_interneuron_paired(df, pyr_col="pv_corr_pyramidal", int_col="pv_corr_interneuron"):
    """
    Paired slope plot for Analysis D: pyramidal-ensemble PV correlation vs.
    interneuron-ensemble PV correlation for the SAME session-pair, each pair
    drawn as two connected points. Makes it directly visible, pair-by-pair,
    whether interneuron PV correlation runs systematically higher / lower /
    more stable than the pyramidal PV correlation for that exact pair --
    independent of, and complementary to, the separate NREM-correlation
    comparison (fig_main_scatter reused with y=pv_corr_interneuron).
    """
    d = df.dropna(subset=[pyr_col, int_col]).reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(5.5, 4.4))
    for _, row in d.iterrows():
        color = color_for_subject(row["subject"])
        ax.plot([0, 1], [row[pyr_col], row[int_col]], color=color, lw=1.3, alpha=0.85,
                marker="o", markersize=6, markeredgecolor="white", markeredgewidth=0.5,
                label=f'{row["subject"]} {row["day_a"]}->{row["day_b"]}')
    ax.set_xlim(-0.3, 1.3)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["Pyramidal\nensemble", "Interneuron\nensemble"])
    ax.set_ylabel("PV correlation")
    ax.set_title(f"Pyramidal vs. interneuron PV correlation, same pairs (n={len(d)})", fontsize=10)
    ax.axhline(0, color="0.7", lw=0.8, ls="--")
    if len(d) <= 8 and len(d) > 0:
        ax.legend(frameon=False, fontsize=6.5, loc="best")
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    return fig


def fig_000061_sleep_architecture(daily_df):
    """
    Descriptive-only figure for DANDI:000061 (NOT used in any PV analysis --
    see notebook). Stacked bars of NREM/REM/Wake duration per recording day,
    grouped by rat, from `processing/behavior/states` (already in seconds;
    see cache_session_000061_descriptive.py for the epochs-vs-states
    time-unit-mismatch investigation).
    """
    daily_df = daily_df.sort_values(["subject", "date"])
    fig, ax = plt.subplots(figsize=(10, 4))
    subjects = daily_df["subject"].unique()
    x = np.arange(len(daily_df))
    state_colors = {"Awake": "#0072B2", "Non-REM": "#D55E00", "REM": "#009E73"}
    bottom = np.zeros(len(daily_df))
    for state in ["Awake", "Non-REM", "REM"]:
        vals = daily_df.get(f"{state}_min", pd.Series(np.zeros(len(daily_df)))).values
        ax.bar(x, vals, bottom=bottom, color=state_colors[state], label=state, width=0.8)
        bottom += vals
    # subject group boundaries
    subj_of_x = daily_df["subject"].values
    boundaries = [i for i in range(1, len(subj_of_x)) if subj_of_x[i] != subj_of_x[i - 1]]
    for b in boundaries:
        ax.axvline(b - 0.5, color="0.7", lw=0.8)
    tick_pos, tick_lab = [], []
    for s in subjects:
        idxs = np.where(subj_of_x == s)[0]
        tick_pos.append(idxs.mean())
        tick_lab.append(s)
    ax.set_xticks(tick_pos); ax.set_xticklabels(tick_lab)
    ax.set_ylabel("Duration (min)")
    ax.set_title("DANDI:000061 -- sleep architecture per recording day (descriptive only, not used for PV)", fontsize=10)
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    return fig
