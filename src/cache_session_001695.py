"""
Cache the essential arrays from one DANDI:001695 (Gonzalez lab) maze session
to a local .npz file, streaming only what is needed via remfile/h5py.

Produces a cache SCHEMA-COMPATIBLE with the DANDI:000213 cache
(cache_session.py) wherever possible -- same field names for position
(pos_t/pos_xy), spikes (spikes_concat/spikes_bounds/n_units), sleep-state
epochs (ss_label/ss_start/ss_stop, using the SAME label vocabulary
{"Awake","Non-REM","REM"} as 000213), and PRE/maze/POST epoch boundaries
(ep_labels/ep_start/ep_stop, using the SAME labels {"before trials",
"during trials","after trials"} as 000213). This lets `analysis_lib.py`
and `build_pairs.py`'s existing `epoch_window`, `state_duration_in_window`,
`is_usable`, `occupancy_and_spike_maps`, `ensemble_tuning_curve`,
`pv_correlation`, and `rv_coefficient` be reused UNMODIFIED on 001695 data
(only the pair-assembly logic -- which environment-pairs are legal, given
this dataset's lack of an explicit epoch table -- needed a new function,
`build_pair_table_001695()` in build_pairs.py).

THREE DOCUMENTED, NON-OBVIOUS DECISIONS MADE HERE (all verified against the
actual data before being adopted, not assumed):

1. CELL SELECTION: this dataset's `units` table has a real waveform-based
   `cell_type` column ({'Pyramidal Cell','Narrow Interneuron',
   'Wide Interneuron','Unknown'}). Only `cell_type == 'Pyramidal Cell'`
   units are kept as place-cell candidates for the PRIMARY (pyramidal)
   ensemble used throughout Analyses B/C. This is a methodological
   IMPROVEMENT over the 000213 pipeline, which has no cell-type column and
   instead falls back on a firing-rate heuristic
   (`min_rate_hz`/`max_rate_hz` in `ensemble_tuning_curve`) to approximate
   "putative pyramidal cell" -- a coarser, unvalidated proxy. Here we use
   the dataset's own ground-truth classification directly.

   ADDED FOR ANALYSIS D: the interneuron set (`cell_type` in
   {'Narrow Interneuron','Wide Interneuron'}) is ALSO cached, in parallel
   fields `n_units_int`/`spikes_concat_int`/`spikes_bounds_int` (same
   concatenated-spike-times layout as the pyramidal `n_units`/
   `spikes_concat`/`spikes_bounds` fields), plus the raw per-class counts
   `n_narrow_interneuron`/`n_wide_interneuron`. This does NOT alter the
   pyramidal fields or any Analysis B/C number -- it is purely additive,
   caching a second, independent unit set from the same session so the
   interneuron ensemble can be run through the exact same
   `occupancy_and_spike_maps`/`ensemble_tuning_curve`/`pv_correlation`
   functions used for pyramidal cells (see `build_pairs.py`,
   `session_pv_inputs_interneuron`).

2. EPOCH INFERENCE (no explicit PRE/maze/POST table exists in this
   dataset -- verified: `intervals` only ever contains `SleepStates`,
   never an `epochs` table). Investigation (see session_list_001695.py
   docstring and notebook Analysis B) found that `processing/behavior/
   AnimalPosition/Position` timestamps span only a SUBSET of the full
   recorded session (e.g. one checked session: position t in
   [4357, 6262] s while sleep-state rows span the full [1, 11858] s) --
   i.e. position tracking is only active DURING the maze bout, not during
   PRE/POST homecage rest. This is used as the epoch-boundary proxy:
     - "during trials" (maze)  = [pos_t.min(), pos_t.max()]
     - "before trials" (PRE)   = [min sleep-state start_time, pos_t.min())
     - "after trials"  (POST)  = (pos_t.max(), max sleep-state stop_time]
   LIMITATION (disclosed prominently in the notebook, not hidden): this
   assumes position tracking is OFF for the entire PRE/POST window and ON
   for the entire maze window, with no gaps. It is a proxy, not an
   experimenter-labelled ground truth epoch table like 000213 has. A
   session where this assumption looks violated (e.g. position tracking
   covering >90% of the full session, leaving no meaningful PRE/POST
   window) is flagged via `epoch_inference_suspect` rather than silently
   accepted.

3. SLEEP-STATE TABLE FILTERING: `intervals/SleepStates` in this dataset
   contains FOUR distinct `state` values: NREM, REM, WAKE, and Ripple.
   Verified directly: 'Ripple' rows have a mean duration of ~26
   milliseconds (vs. tens-to-hundreds of seconds for NREM/REM/WAKE) --
   these are individual sharp-wave-ripple event timestamps, NOT a
   behavioral state occupying continuous time, and are DROPPED before any
   "between-session NREM" duration is computed (they would otherwise
   silently inflate row counts without representing sleep time). NREM/REM/
   WAKE are relabelled Non-REM/REM/Awake respectively to match the 000213
   vocabulary. Ripple event count/rate is still cached (not discarded) for
   independent QC/diagnostic use, just not fed into `state_duration_in_window`.

Usage: python3 cache_session_001695.py <path_in_dandiset> <outfile.npz>
"""
import sys
import numpy as np
import remfile, h5py, requests


def get_download_url(dandiset_id, path, version="draft"):
    r = requests.get(
        f"https://api.dandiarchive.org/api/dandisets/{dandiset_id}/versions/{version}/assets/",
        params={"path": path}, timeout=60,
    )
    r.raise_for_status()
    d = r.json()
    if not d["results"]:
        raise FileNotFoundError(path)
    asset_id = d["results"][0]["asset_id"]
    r2 = requests.get(
        f"https://api.dandiarchive.org/api/dandisets/{dandiset_id}/versions/{version}/assets/{asset_id}/",
        timeout=60,
    )
    r2.raise_for_status()
    urls = r2.json()["contentUrl"]
    return urls[1] if len(urls) > 1 else urls[0]


def _decode(x):
    return x.decode() if isinstance(x, bytes) else x


STATE_RELABEL = {"WAKE": "Awake", "NREM": "Non-REM", "REM": "REM"}  # 'Ripple' dropped, not relabelled


def cache_session(dandiset_id, path, outfile):
    url = get_download_url(dandiset_id, path)
    rf = remfile.File(url)
    f = h5py.File(rf, "r")

    session_start_time = _decode(f["session_start_time"][()]) if "session_start_time" in f else ""
    session_description = _decode(f["session_description"][()]) if "session_description" in f else ""
    subject_id = _decode(f["general/subject/subject_id"][()]) if "general/subject/subject_id" in f else ""

    # --- position (required for this analysis; caller should skip sessions without it) ---
    pos_path = "processing/behavior/AnimalPosition/Position"
    has_position = pos_path in f
    if has_position:
        pos = f[pos_path]
        pos_t_raw = pos["timestamps"][:].astype(np.float64)
        pos_xy_raw = pos["data"][:, :2].astype(np.float64)
        finite = np.isfinite(pos_t_raw) & np.isfinite(pos_xy_raw).all(axis=1)
        pos_t = pos_t_raw[finite]
        pos_xy = pos_xy_raw[finite]
    else:
        pos_t = np.array([])
        pos_xy = np.array([]).reshape(0, 2)

    # --- speed (cached for QC/diagnostics; not required by the core pipeline
    #     since epoch inference uses the position-coverage window directly) ---
    spd_path = "processing/behavior/Speed"
    has_speed = spd_path in f
    if has_speed:
        spd = f[spd_path]
        spd_t = spd["timestamps"][:].astype(np.float64)
        spd_v = spd["data"][:].astype(np.float64)
    else:
        spd_t = np.array([])
        spd_v = np.array([])

    # --- sleep states (NREM/REM/WAKE kept + relabelled; Ripple dropped) ---
    ss_path = "intervals/SleepStates"
    has_sleep_states = ss_path in f
    n_ripple_events = 0
    if has_sleep_states:
        ss = f[ss_path]
        raw_label = np.array([_decode(x) for x in ss["state"][:]])
        raw_start = ss["start_time"][:].astype(float)
        raw_stop = ss["stop_time"][:].astype(float)
        n_ripple_events = int((raw_label == "Ripple").sum())
        keep = np.isin(raw_label, list(STATE_RELABEL.keys()))
        ss_label = np.array([STATE_RELABEL[x] for x in raw_label[keep]])
        ss_start = raw_start[keep]
        ss_stop = raw_stop[keep]
    else:
        ss_label, ss_start, ss_stop = np.array([]), np.array([]), np.array([])

    # --- epoch inference from position coverage (see module docstring #2) ---
    has_epochs = False
    ep_labels, ep_start, ep_stop = [], [], []
    epoch_inference_suspect = True
    if has_position and len(pos_t) > 10 and has_sleep_states and len(ss_start) > 0:
        maze_t0, maze_t1 = float(pos_t.min()), float(pos_t.max())
        sess_t0, sess_t1 = float(ss_start.min()), float(ss_stop.max())
        pre_dur = max(0.0, maze_t0 - sess_t0)
        post_dur = max(0.0, sess_t1 - maze_t1)
        maze_dur = maze_t1 - maze_t0
        total_dur = max(sess_t1 - sess_t0, 1e-9)
        # Flag as suspect if the inferred maze window swallows almost the
        # whole session (leaves < 2% of total time as PRE+POST combined) --
        # in that case PRE/POST rest cannot be meaningfully separated from
        # the maze bout, and downstream between-session sleep numbers for
        # this session should be treated with extra caution.
        epoch_inference_suspect = (pre_dur + post_dur) < 0.02 * total_dur
        ep_labels = ["before trials", "during trials", "after trials"]
        ep_start = [sess_t0, maze_t0, maze_t1]
        ep_stop = [maze_t0, maze_t1, sess_t1]
        has_epochs = True

    # --- units: cell_type == 'Pyramidal Cell' only (see module docstring #1) ---
    # ALSO cache the interneuron set (Narrow Interneuron + Wide Interneuron)
    # in parallel fields (n_units_int/spikes_concat_int/spikes_bounds_int),
    # for the Analysis D interneuron-vs-pyramidal comparison (see notebook).
    # This does NOT change which units feed the pyramidal-ensemble pipeline
    # used by Analyses B/C -- n_units/spikes_concat/spikes_bounds below are
    # unchanged pyramidal-only fields, exactly as before.
    n_units = 0
    spikes_concat = np.array([], dtype=np.float64)
    spikes_bounds = np.array([0], dtype=np.int64)
    n_units_int = 0
    spikes_concat_int = np.array([], dtype=np.float64)
    spikes_bounds_int = np.array([0], dtype=np.int64)
    n_narrow_interneuron = 0
    n_wide_interneuron = 0
    has_cell_type_col = "units" in f and "cell_type" in f["units"]
    if has_cell_type_col:
        u = f["units"]
        ctype = np.array([_decode(x) for x in u["cell_type"][:]])
        sp_idx_all = u["spike_times_index"][:]
        sp_data_all = u["spike_times"]

        def _gather(idx_arr):
            spike_list = []
            for i in idx_arr:
                lo = sp_idx_all[i - 1] if i > 0 else 0
                hi = sp_idx_all[i]
                spike_list.append(np.asarray(sp_data_all[lo:hi], dtype=np.float64))
            n = len(spike_list)
            if n > 0:
                concat = np.concatenate(spike_list)
                bounds = np.concatenate([[0], np.cumsum([len(s) for s in spike_list])])
            else:
                concat = np.array([], dtype=np.float64)
                bounds = np.array([0], dtype=np.int64)
            return n, concat, bounds

        pyr_idx = np.where(ctype == "Pyramidal Cell")[0]
        n_units, spikes_concat, spikes_bounds = _gather(pyr_idx)

        int_idx = np.where(np.isin(ctype, ["Narrow Interneuron", "Wide Interneuron"]))[0]
        n_units_int, spikes_concat_int, spikes_bounds_int = _gather(int_idx)
        n_narrow_interneuron = int((ctype == "Narrow Interneuron").sum())
        n_wide_interneuron = int((ctype == "Wide Interneuron").sum())

    # --- footprint (for same-environment clustering; see session_list_001695.py) ---
    if has_position and len(pos_xy) > 0:
        x_span = float(np.nanmax(pos_xy[:, 0]) - np.nanmin(pos_xy[:, 0]))
        y_span = float(np.nanmax(pos_xy[:, 1]) - np.nanmin(pos_xy[:, 1]))
    else:
        x_span, y_span = np.nan, np.nan

    f.close()

    np.savez_compressed(
        outfile,
        session_start_time=session_start_time,
        session_description=session_description,
        subject_id=subject_id,
        path=path,
        has_position=has_position,
        has_speed=has_speed,
        has_sleep_states=has_sleep_states,
        has_epochs=has_epochs,
        epoch_inference_suspect=epoch_inference_suspect,
        has_cell_type_col=has_cell_type_col,
        ep_labels=np.array(ep_labels), ep_start=np.array(ep_start), ep_stop=np.array(ep_stop),
        ss_label=ss_label, ss_start=ss_start, ss_stop=ss_stop,
        n_ripple_events=n_ripple_events,
        n_units=n_units,
        spikes_concat=spikes_concat, spikes_bounds=spikes_bounds,
        n_units_int=n_units_int,
        spikes_concat_int=spikes_concat_int, spikes_bounds_int=spikes_bounds_int,
        n_narrow_interneuron=n_narrow_interneuron, n_wide_interneuron=n_wide_interneuron,
        pos_t=pos_t, pos_xy=pos_xy,
        spd_t=spd_t, spd_v=spd_v,
        x_span=x_span, y_span=y_span,
    )
    usable = has_epochs and has_sleep_states and has_position and n_units >= 5
    print(f"OK {outfile} subj={subject_id} n_pyr_units={n_units} n_int_units={n_units_int} "
          f"(narrow={n_narrow_interneuron} wide={n_wide_interneuron}) n_pos={len(pos_t)} "
          f"has_epochs={has_epochs} suspect={epoch_inference_suspect} usable={usable}", flush=True)


if __name__ == "__main__":
    dandiset_id = "001695"
    path = sys.argv[1]
    outfile = sys.argv[2]
    try:
        cache_session(dandiset_id, path, outfile)
    except Exception as e:
        print(f"FAIL {path}: {e}", flush=True)
        sys.exit(1)
