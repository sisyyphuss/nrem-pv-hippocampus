"""
Core analysis functions for the NREM-sleep vs. hippocampal population-vector
(PV) correlation pipeline (DANDI:000213, Tingley lab, Buzsaki lab).

All functions operate on the cached per-session .npz arrays produced by
cache_session.py. Kept dependency-light (numpy/scipy only) so it is easy to
audit every computation by eye.
"""
import numpy as np
from scipy.stats import zscore


# ----------------------------------------------------------------------
# Session-level sleep quantification
# ----------------------------------------------------------------------

def epoch_window(sess, name):
    """Return (start, stop) for an epoch label ('before trials'/'during trials'/'after trials')."""
    labels = sess["ep_labels"]
    idx = np.where(labels == name)[0]
    if len(idx) == 0:
        return None
    i = idx[0]
    return float(sess["ep_start"][i]), float(sess["ep_stop"][i])


def state_duration_in_window(sess, state_label, t0, t1):
    """Total seconds spent in `state_label` (e.g. 'Non-REM') within [t0, t1)."""
    lab = sess["ss_label"]
    s0 = sess["ss_start"].astype(float)
    s1 = sess["ss_stop"].astype(float)
    mask = lab == state_label
    lo = np.maximum(s0[mask], t0)
    hi = np.minimum(s1[mask], t1)
    dur = np.clip(hi - lo, 0, None)
    return float(dur.sum())


def post_trial_sleep(sess, state_label="Non-REM"):
    """Duration of `state_label` during the 'after trials' epoch (this day's POST-behavior rest)."""
    win = epoch_window(sess, "after trials")
    if win is None:
        return np.nan
    return state_duration_in_window(sess, state_label, *win)


def pre_trial_sleep(sess, state_label="Non-REM"):
    """Duration of `state_label` during the 'before trials' epoch (this day's PRE-behavior rest)."""
    win = epoch_window(sess, "before trials")
    if win is None:
        return np.nan
    return state_duration_in_window(sess, state_label, *win)


def between_session_state(sess_day_n, sess_day_n1, state_label="Non-REM"):
    """
    Best available proxy for time spent in `state_label` *between* two
    recording sessions: POST-trial rest of day N + PRE-trial rest of day N+1.
    NOTE: this necessarily misses any untracked overnight sleep that falls
    outside the recorded epochs (see notebook Limitations section).
    """
    post_n = post_trial_sleep(sess_day_n, state_label)
    pre_n1 = pre_trial_sleep(sess_day_n1, state_label)
    return post_n + pre_n1


# ----------------------------------------------------------------------
# Place fields / rate maps
# ----------------------------------------------------------------------

def get_unit_spikes(sess, unit_idx):
    b = sess["spikes_bounds"]
    return sess["spikes_concat"][b[unit_idx]:b[unit_idx + 1]]


def compute_speed(pos_t, pos_xy):
    dt = np.diff(pos_t)
    dt[dt <= 0] = np.nan
    dxy = np.diff(pos_xy, axis=0)
    dist = np.sqrt((dxy ** 2).sum(axis=1))
    speed = dist / dt
    speed = np.concatenate([[speed[0]], speed])
    return speed


def occupancy_and_spike_maps(sess, t0, t1, bin_edges_x, bin_edges_y,
                              speed_thresh=None, min_occ_sec=0.15):
    """
    Build 2D occupancy map (seconds/bin) restricted to [t0,t1) and, for every
    hippocampal unit, its spike-count map -> returns rate maps (Hz),
    occupancy mask (bins with enough sampling), and n_units.
    """
    pos_t = sess["pos_t"]
    pos_xy = sess["pos_xy"]
    m = (pos_t >= t0) & (pos_t < t1)
    pt, pxy = pos_t[m], pos_xy[m]
    if len(pt) < 10:
        return None

    if speed_thresh is not None:
        spd = compute_speed(pt, pxy)
        run_mask = spd >= speed_thresh
        pt, pxy = pt[run_mask], pxy[run_mask]
        if len(pt) < 10:
            return None

    # per-sample dt for occupancy (irregular sampling)
    dt = np.diff(pt, append=pt[-1] + np.median(np.diff(pt)))
    dt = np.clip(dt, 0, 0.5)  # guard against big position gaps inflating occupancy

    occ, _, _ = np.histogram2d(pxy[:, 0], pxy[:, 1],
                                bins=[bin_edges_x, bin_edges_y], weights=dt)
    occ_mask = occ >= min_occ_sec

    n_units = int(sess["n_units"])
    rate_maps = np.full((n_units, occ.shape[0], occ.shape[1]), np.nan)
    for u in range(n_units):
        sp = get_unit_spikes(sess, u)
        sp = sp[(sp >= t0) & (sp < t1)]
        if speed_thresh is not None and len(sp) > 0:
            # keep only spikes that fall within a position sample already speed-filtered
            # (nearest-neighbour lookup on filtered position timestamps)
            idxs = np.searchsorted(pt, sp)
            idxs = np.clip(idxs, 0, len(pt) - 1)
            keep = np.abs(pt[idxs] - sp) < 0.25
            sp = sp[keep]
        if len(sp) == 0:
            spike_xy = np.empty((0, 2))
        else:
            idxs = np.searchsorted(pt, sp)
            idxs = np.clip(idxs, 0, len(pt) - 1)
            spike_xy = pxy[idxs]
        sc, _, _ = np.histogram2d(spike_xy[:, 0], spike_xy[:, 1],
                                   bins=[bin_edges_x, bin_edges_y])
        with np.errstate(divide="ignore", invalid="ignore"):
            rm = sc / occ
        rm[~occ_mask] = np.nan
        rate_maps[u] = rm

    return {"rate_maps": rate_maps, "occupancy": occ, "occ_mask": occ_mask}


def ensemble_tuning_curve(rate_maps, occ_mask, min_rate_hz=0.1, max_rate_hz=None):
    """
    Collapse a [n_units, nx, ny] rate-map stack into a single ensemble tuning
    curve: per-unit z-score across valid bins, then average across units.
    Units with near-zero overall firing are excluded (silent/too-sparse cells
    give an undefined z-score and would inject noise).
    Returns (tuning_vector[valid_bins], n_units_used).
    """
    nU = rate_maps.shape[0]
    flat_mask = occ_mask.ravel()
    curves = []
    for u in range(nU):
        rm = rate_maps[u].ravel()[flat_mask]
        if np.nanmean(rm) < min_rate_hz:
            continue
        if max_rate_hz is not None and np.nanmean(rm) > max_rate_hz:
            continue
        rm = np.nan_to_num(rm, nan=0.0)
        if rm.std() == 0:
            continue
        curves.append(zscore(rm))
    if len(curves) == 0:
        return None, 0
    curves = np.array(curves)
    return curves.mean(axis=0), len(curves)


def pv_correlation(vec1, vec2):
    """Pearson correlation between two ensemble tuning-curve vectors (same bin order)."""
    if vec1 is None or vec2 is None:
        return np.nan
    ok = np.isfinite(vec1) & np.isfinite(vec2)
    if ok.sum() < 5:
        return np.nan
    r = np.corrcoef(vec1[ok], vec2[ok])[0, 1]
    return r


def rv_coefficient(mat1, mat2, occ_mask1, occ_mask2):
    """
    RV coefficient: a multivariate similarity between two [bins x cells]
    matrices that does NOT require matched cell identity or matched cell
    count between sessions (unlike a literal per-cell PV correlation).
    Used as an independent-method robustness check on pv_correlation().
    mat1/mat2: [n_units, nx, ny] rate maps (as from occupancy_and_spike_maps).
    """
    common_mask = occ_mask1.ravel() & occ_mask2.ravel()
    if common_mask.sum() < 5:
        return np.nan

    def prep(mat, mask):
        nU = mat.shape[0]
        X = np.stack([mat[u].ravel()[mask] for u in range(nU)], axis=1)  # bins x cells
        X = np.nan_to_num(X, nan=0.0)
        keep = X.std(axis=0) > 0
        X = X[:, keep]
        if X.shape[1] == 0:
            return None
        X = (X - X.mean(axis=0)) / X.std(axis=0)
        return X

    X1 = prep(mat1, common_mask)
    X2 = prep(mat2, common_mask)
    if X1 is None or X2 is None:
        return np.nan
    S1 = X1 @ X1.T  # bins x bins similarity
    S2 = X2 @ X2.T
    num = np.trace(S1 @ S2)
    den = np.sqrt(np.trace(S1 @ S1) * np.trace(S2 @ S2))
    if den == 0:
        return np.nan
    return float(num / den)
