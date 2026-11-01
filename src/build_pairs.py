"""
Build the session-pair dataset: for every same-subject, same-task pair of
sessions (day N, day N+1-in-task-sequence), compute:
  - the PV correlation (ensemble tuning-curve method + RV-coefficient method)
  - the amount of NREM / REM / Wake recorded between the two sessions
  - the calendar-day gap and an explicit flag for whether any intervening
    recorded (but excluded) days exist, and how much of their sleep is
    folded into the "between-session" sleep total
"""
import glob, os
import numpy as np
import pandas as pd
from datetime import datetime

from analysis_lib import (
    epoch_window, state_duration_in_window, occupancy_and_spike_maps,
    ensemble_tuning_curve, pv_correlation, rv_coefficient,
)
from session_list import SESSIONS


def load_all_cached():
    """Load every cached .npz (including off-task / zero-unit ones) keyed by subject+path."""
    out = {}
    for subj, sessions in SESSIONS.items():
        out[subj] = []
        for idx, path in enumerate(sessions):
            f = f"data/{subj}_{idx}.npz"
            if not os.path.exists(f):
                continue
            d = np.load(f, allow_pickle=True)
            rec = {k: d[k] for k in d.files}
            rec["date"] = datetime.fromisoformat(str(rec["session_start_time"])[:19])
            out[subj].append(rec)
        out[subj].sort(key=lambda r: r["date"])
    return out


def total_day_state_duration(sess, state_label="Non-REM"):
    if sess["ss_label"].size == 0:
        return np.nan
    lab = sess["ss_label"]
    dur = (sess["ss_stop"].astype(float) - sess["ss_start"].astype(float))[lab == state_label]
    return float(dur.sum())


def between_session_sleep(all_sessions_for_subject, i_from, i_to, state_label="Non-REM"):
    """
    NREM (or REM/Wake) recorded between session index i_from and i_to in the
    subject's FULL chronological session list (which may include off-task /
    zero-unit days that were still recorded and still contain valid sleep
    scoring). Sums: POST-trial rest of day i_from + full-day sleep of any
    strictly intervening recorded days + PRE-trial rest of day i_to.
    """
    s_from = all_sessions_for_subject[i_from]
    s_to = all_sessions_for_subject[i_to]

    total = 0.0
    win = epoch_window(s_from, "after trials")
    if win is not None and s_from["ss_label"].size:
        total += state_duration_in_window(s_from, state_label, *win)

    for j in range(i_from + 1, i_to):
        total += total_day_state_duration(all_sessions_for_subject[j], state_label)

    win2 = epoch_window(s_to, "before trials")
    if win2 is not None and s_to["ss_label"].size:
        total += state_duration_in_window(s_to, state_label, *win2)

    return total


def is_usable(sess, min_units=5):
    return (
        bool(sess["has_epochs"]) and bool(sess["has_sleep_states"]) and
        bool(sess["has_position"]) and int(sess["n_units"]) >= min_units
    )


def build_shared_bins(sess_a, sess_b, n_bins=12):
    """Common 2D bin edges spanning both sessions' explored positions."""
    xy = np.concatenate([sess_a["pos_xy"], sess_b["pos_xy"]], axis=0)
    xy = xy[np.isfinite(xy).all(axis=1)]
    x_edges = np.linspace(np.percentile(xy[:, 0], 0.5), np.percentile(xy[:, 0], 99.5), n_bins + 1)
    y_edges = np.linspace(np.percentile(xy[:, 1], 0.5), np.percentile(xy[:, 1], 99.5), n_bins + 1)
    return x_edges, y_edges


def session_pv_inputs(sess, x_edges, y_edges, speed_thresh=None):
    win = epoch_window(sess, "during trials")
    if win is None:
        return None
    t0, t1 = win
    return occupancy_and_spike_maps(sess, t0, t1, x_edges, y_edges, speed_thresh=speed_thresh)


def build_pair_table(n_bins=12, speed_thresh=None, min_units=5, max_day_gap=None):
    all_sessions = load_all_cached()
    rows = []
    for subj, sessions in all_sessions.items():
        n = len(sessions)
        # restrict to the task used on the FIRST usable day of this subject,
        # to define "the familiar environment" for that animal
        usable_idx = [i for i in range(n) if is_usable(sessions[i], min_units)]
        if not usable_idx:
            continue
        task_counts = {}
        for i in usable_idx:
            t = str(sessions[i]["task_group"])
            task_counts[t] = task_counts.get(t, 0) + 1
        main_task = max(task_counts, key=task_counts.get)
        same_task_idx = [i for i in usable_idx if str(sessions[i]["task_group"]) == main_task]

        for a, b in zip(same_task_idx[:-1], same_task_idx[1:]):
            sess_a, sess_b = sessions[a], sessions[b]
            day_gap = (sess_b["date"] - sess_a["date"]).days
            if max_day_gap is not None and day_gap > max_day_gap:
                continue

            x_edges, y_edges = build_shared_bins(sess_a, sess_b, n_bins=n_bins)
            in_a = session_pv_inputs(sess_a, x_edges, y_edges, speed_thresh)
            in_b = session_pv_inputs(sess_b, x_edges, y_edges, speed_thresh)
            if in_a is None or in_b is None:
                continue

            # Use a SHARED valid-bin mask (bins with enough occupancy in BOTH
            # sessions) so the two ensemble tuning curves are the same length
            # and directly comparable bin-for-bin. Using each session's own
            # mask independently (as an earlier version of this function did)
            # produces mismatched vector lengths whenever the two sessions'
            # explored/undersampled bins differ -- caught by execution, not
            # left silent.
            common_mask = in_a["occ_mask"] & in_b["occ_mask"]
            vec_a, nu_a = ensemble_tuning_curve(in_a["rate_maps"], common_mask)
            vec_b, nu_b = ensemble_tuning_curve(in_b["rate_maps"], common_mask)
            if vec_a is None or vec_b is None or nu_a < 3 or nu_b < 3:
                continue

            pvr = pv_correlation(vec_a, vec_b)
            rvc = rv_coefficient(in_a["rate_maps"], in_b["rate_maps"], in_a["occ_mask"], in_b["occ_mask"])

            nrem_between = between_session_sleep(sessions, a, b, "Non-REM")
            rem_between = between_session_sleep(sessions, a, b, "REM")
            wake_between = between_session_sleep(sessions, a, b, "Awake")

            rows.append(dict(
                subject=subj, task=main_task,
                day_a=sess_a["date"].date().isoformat(), day_b=sess_b["date"].date().isoformat(),
                day_gap=day_gap,
                n_units_a=int(sess_a["n_units"]), n_units_b=int(sess_b["n_units"]),
                n_units_used_a=nu_a, n_units_used_b=nu_b,
                pv_corr=pvr, pv_corr_rv=rvc,
                nrem_between_min=nrem_between / 60.0,
                rem_between_min=rem_between / 60.0,
                wake_between_min=wake_between / 60.0,
            ))
    return pd.DataFrame(rows)


if __name__ == "__main__":
    df = build_pair_table()
    print(df.to_string())
    df.to_csv("session_pairs.csv", index=False)
    print(f"\n{len(df)} session-pairs written to session_pairs.csv")


# ============================================================================
# EXTENSION FOR DANDI:001695 (Analysis B / C).
#
# Added as NEW functions rather than modifying build_pair_table() above,
# per the project convention of not rewriting existing, already-validated
# 000213-specific code. The 001695 cache (cache_session_001695.py) was
# deliberately built to reuse the SAME field names and label vocabulary as
# the 000213 cache (ep_labels/ep_start/ep_stop with "before/during/after
# trials"; ss_label with "Awake"/"Non-REM"/"REM"), so `epoch_window`,
# `state_duration_in_window`, `is_usable`, `build_shared_bins`,
# `session_pv_inputs`, `ensemble_tuning_curve`, `pv_correlation`, and
# `rv_coefficient` are reused UNMODIFIED below -- only the pair-selection
# logic (WHICH sessions are legal to pair) differs, because 001695 has no
# task-group label to key off of; it uses the aspect-ratio same-environment
# clustering documented in session_list_001695.py instead.
# ============================================================================

from session_list_001695 import SESSIONS_001695, ASPECT_RATIO_THRESHOLD


def load_all_cached_001695():
    """Load every cached 001695 .npz, keyed by subject, chronologically sorted."""
    out = {}
    for subj, sessions in SESSIONS_001695.items():
        out[subj] = []
        for idx, path in enumerate(sessions):
            f = f"data_001695/{subj}_{idx}.npz"
            if not os.path.exists(f):
                continue
            d = np.load(f, allow_pickle=True)
            rec = {k: d[k] for k in d.files}
            rec["date"] = datetime.fromisoformat(str(rec["session_start_time"])[:19])
            out[subj].append(rec)
        out[subj].sort(key=lambda r: r["date"])
    return out


def maze_config_label(sess, threshold=ASPECT_RATIO_THRESHOLD):
    """
    'narrow' vs 'wide' apparatus-configuration label from the position
    footprint aspect ratio (y_span / x_span), using the data-driven
    threshold derived in session_list_001695.py. Returns None if the
    footprint could not be computed (e.g. no position data).
    """
    x_span, y_span = float(sess["x_span"]), float(sess["y_span"])
    if not (np.isfinite(x_span) and np.isfinite(y_span)) or x_span <= 0:
        return None
    return "wide" if (y_span / x_span) >= threshold else "narrow"


def build_pair_table_001695(n_bins=12, speed_thresh=None, min_units=5):
    """
    Same-environment session pairs for DANDI:001695. A pair (day N, day M)
    is included only if BOTH sessions carry the SAME maze_config_label
    (see session_list_001695.py for why: the uniform "Novel Maze Session"
    text label cannot distinguish reused-vs-new environments, so the
    apparatus footprint is used as the actual evidence). Within a same-
    label run, only CONSECUTIVE sessions in that run are paired -- directly
    analogous to build_pair_table()'s same-task consecutive-pairing logic
    for 000213.
    """
    all_sessions = load_all_cached_001695()
    rows = []
    for subj, sessions in all_sessions.items():
        n = len(sessions)
        usable_idx = [i for i in range(n) if is_usable(sessions[i], min_units)]
        if not usable_idx:
            continue

        labels = {i: maze_config_label(sessions[i]) for i in usable_idx}
        # group indices into runs of consecutive (in the usable, chronological
        # sequence) sessions sharing the same maze_config_label
        clusters = {}
        for i in usable_idx:
            clusters.setdefault(labels[i], []).append(i)

        for cluster_label, idx_list in clusters.items():
            if cluster_label is None or len(idx_list) < 2:
                continue
            idx_list = sorted(idx_list)
            for a, b in zip(idx_list[:-1], idx_list[1:]):
                sess_a, sess_b = sessions[a], sessions[b]
                day_gap = (sess_b["date"] - sess_a["date"]).days

                x_edges, y_edges = build_shared_bins(sess_a, sess_b, n_bins=n_bins)
                in_a = session_pv_inputs(sess_a, x_edges, y_edges, speed_thresh)
                in_b = session_pv_inputs(sess_b, x_edges, y_edges, speed_thresh)
                if in_a is None or in_b is None:
                    continue

                common_mask = in_a["occ_mask"] & in_b["occ_mask"]
                vec_a, nu_a = ensemble_tuning_curve(in_a["rate_maps"], common_mask)
                vec_b, nu_b = ensemble_tuning_curve(in_b["rate_maps"], common_mask)
                if vec_a is None or vec_b is None or nu_a < 3 or nu_b < 3:
                    continue

                pvr = pv_correlation(vec_a, vec_b)
                rvc = rv_coefficient(in_a["rate_maps"], in_b["rate_maps"], in_a["occ_mask"], in_b["occ_mask"])

                nrem_between = between_session_sleep(sessions, a, b, "Non-REM")
                rem_between = between_session_sleep(sessions, a, b, "REM")
                wake_between = between_session_sleep(sessions, a, b, "Awake")

                rows.append(dict(
                    dataset="001695", species="mouse",
                    subject=subj, task=f"maze_config:{cluster_label}",
                    day_a=sess_a["date"].date().isoformat(), day_b=sess_b["date"].date().isoformat(),
                    day_gap=day_gap,
                    n_units_a=int(sess_a["n_units"]), n_units_b=int(sess_b["n_units"]),
                    n_units_used_a=nu_a, n_units_used_b=nu_b,
                    pv_corr=pvr, pv_corr_rv=rvc,
                    nrem_between_min=nrem_between / 60.0,
                    rem_between_min=rem_between / 60.0,
                    wake_between_min=wake_between / 60.0,
                    epoch_inference_suspect=bool(sess_a["epoch_inference_suspect"]) or bool(sess_b["epoch_inference_suspect"]),
                ))
    return pd.DataFrame(rows)


# ============================================================================
# EXTENSION FOR ANALYSIS D: interneuron vs. pyramidal sleep-dependent drift
# dissociation (DANDI:001695 only -- 000213's units table has no cell_type
# column, so it cannot do this comparison; see cache_session_001695.py).
#
# Added as NEW functions, per the same project convention as the 001695
# extension above: `occupancy_and_spike_maps`, `ensemble_tuning_curve`, and
# `pv_correlation` are reused COMPLETELY UNMODIFIED. The only new piece of
# logic is `session_pv_inputs_interneuron`, which builds a shallow "shadow"
# session dict pointing the pyramidal-named fields
# (n_units/spikes_concat/spikes_bounds) that occupancy_and_spike_maps()
# expects at the INTERNEURON fields (n_units_int/spikes_concat_int/
# spikes_bounds_int) cached by cache_session_001695.py, so no analysis
# function needs to know "interneuron" exists as a concept at all.
# ============================================================================

def session_pv_inputs_interneuron(sess, x_edges, y_edges, speed_thresh=None):
    """
    Same as session_pv_inputs(), but built on the INTERNEURON unit set
    (cell_type in {'Narrow Interneuron','Wide Interneuron'}) instead of the
    pyramidal set, by swapping in the n_units_int/spikes_concat_int/
    spikes_bounds_int fields under the names occupancy_and_spike_maps()
    expects. Returns None if this session's cache predates Analysis D (no
    n_units_int field) or has no 'during trials' epoch window.
    """
    if "n_units_int" not in sess:
        return None
    win = epoch_window(sess, "during trials")
    if win is None:
        return None
    t0, t1 = win
    shadow = dict(sess)
    shadow["n_units"] = sess["n_units_int"]
    shadow["spikes_concat"] = sess["spikes_concat_int"]
    shadow["spikes_bounds"] = sess["spikes_bounds_int"]
    return occupancy_and_spike_maps(shadow, t0, t1, x_edges, y_edges, speed_thresh=speed_thresh)


def build_pair_table_001695_interneuron(pyramidal_df, n_bins=12, speed_thresh=None):
    """
    For every session-pair ALREADY present in Analysis B's pyramidal pair
    table (`pyramidal_df`, the output of build_pair_table_001695()),
    recompute the ensemble tuning curve / PV correlation for the exact same
    two sessions using ONLY the interneuron unit set. This guarantees Analysis
    D uses the identical set of session-pairs as Analysis B (same subjects,
    same same-environment clustering, same day_a/day_b) -- only the unit
    class fed into ensemble_tuning_curve() differs.

    Does NOT lower ensemble_tuning_curve()'s own minimum-unit floor: a pair
    is only marked `usable_interneuron=True` if at least 3 interneurons
    survive the tuning-curve's internal firing-rate/variance filters on
    BOTH sessions (identical threshold to the pyramidal pipeline's nu>=3
    check in build_pair_table_001695()). Pairs that fail this are kept in
    the returned table (not silently dropped) with pv_corr_interneuron=NaN
    and usable_interneuron=False, so the true attempted-vs-usable
    denominator stays visible, matching this project's honest-reporting
    standard.

    Raw (pre-filter) interneuron counts per session (n_int_units_a/b) are
    reported for every row regardless of usability, for direct inspection
    of whether low counts (rather than the tuning-curve's internal
    firing-rate filter) are the limiting factor.
    """
    all_sessions = load_all_cached_001695()
    rows = []
    for _, r in pyramidal_df.iterrows():
        subj = r["subject"]
        sessions = all_sessions[subj]
        sess_a = next(s for s in sessions if s["date"].date().isoformat() == r["day_a"])
        sess_b = next(s for s in sessions if s["date"].date().isoformat() == r["day_b"])

        n_int_a = int(sess_a["n_units_int"]) if "n_units_int" in sess_a else 0
        n_int_b = int(sess_b["n_units_int"]) if "n_units_int" in sess_b else 0

        x_edges, y_edges = build_shared_bins(sess_a, sess_b, n_bins=n_bins)
        in_a = session_pv_inputs_interneuron(sess_a, x_edges, y_edges, speed_thresh)
        in_b = session_pv_inputs_interneuron(sess_b, x_edges, y_edges, speed_thresh)

        pv_int, nu_a_int, nu_b_int = np.nan, 0, 0
        usable_int = False
        if in_a is not None and in_b is not None:
            common_mask = in_a["occ_mask"] & in_b["occ_mask"]
            vec_a, nu_a_int = ensemble_tuning_curve(in_a["rate_maps"], common_mask)
            vec_b, nu_b_int = ensemble_tuning_curve(in_b["rate_maps"], common_mask)
            if vec_a is not None and vec_b is not None and nu_a_int >= 3 and nu_b_int >= 3:
                pv_int = pv_correlation(vec_a, vec_b)
                usable_int = True

        rows.append(dict(
            subject=subj, task=r["task"],
            day_a=r["day_a"], day_b=r["day_b"], day_gap=r["day_gap"],
            n_int_units_a=n_int_a, n_int_units_b=n_int_b,
            n_int_units_used_a=nu_a_int, n_int_units_used_b=nu_b_int,
            pv_corr_interneuron=pv_int, usable_interneuron=usable_int,
            pv_corr_pyramidal=r["pv_corr"],
            nrem_between_min=r["nrem_between_min"],
            rem_between_min=r["rem_between_min"],
            wake_between_min=r["wake_between_min"],
        ))
    return pd.DataFrame(rows)
