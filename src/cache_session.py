"""
Cache the essential arrays from one Tingley-lab (DANDI:000213) NWB session
to a local .npz file, streaming only what is needed via remfile/h5py
(no full-file download).

Robust to the structural variability discovered during QC: task/apparatus
group name under processing/behavior differs across subjects and even
across days for the same subject (e.g. 'central_alternation_task' vs
'wheel_alternation_task' vs 'linear_task'), and some sessions lack
sleep_states, epochs, or a units 'location' column entirely. All such
gaps are recorded explicitly in the cache (as empty arrays / flags) rather
than silently guessed at, so the notebook's QC step can filter on them.

Usage: python3 cache_session.py <dandiset_path_in_archive> <outfile.npz>
"""
import sys, os, json
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


def cache_session(dandiset_id, path, outfile):
    url = get_download_url(dandiset_id, path)
    rf = remfile.File(url)
    f = h5py.File(rf, "r")

    session_start_time = f.attrs.get(
        "session_start_time", f["session_start_time"][()] if "session_start_time" in f else b""
    )
    if isinstance(session_start_time, bytes):
        session_start_time = session_start_time.decode()
    subject_id = f["general/subject/subject_id"][()].decode()

    # --- task/apparatus group (name varies across sessions -- discover, don't hardcode) ---
    beh_keys = list(f["processing/behavior"].keys()) if "processing/behavior" in f else []
    task_candidates = [k for k in beh_keys if k != "allocentric_frame_tracking"]
    task_group = None
    for k in task_candidates:
        if "position" in f[f"processing/behavior/{k}"]:
            task_group = k
            break

    # --- epochs (before/during/after trials), if present ---
    has_epochs = "intervals/epochs" in f
    ep_labels, ep_start, ep_stop = [], np.array([]), np.array([])
    if has_epochs:
        ep = f["intervals/epochs"]
        tags_flat = ep["tags"][:]
        tags_idx = ep["tags_index"][:]
        ep_start = ep["start_time"][:]
        ep_stop = ep["stop_time"][:]
        prev = 0
        for i in range(len(ep_start)):
            tg = tags_flat[prev:tags_idx[i]]
            prev = tags_idx[i]
            ep_labels.append(tg[0].decode() if len(tg) else "")

    # --- sleep states, if present ---
    has_sleep_states = "processing/ecephys/sleep_states" in f
    if has_sleep_states:
        ss = f["processing/ecephys/sleep_states"]
        ss_label = np.array([x.decode() for x in ss["label"][:]])
        ss_start = ss["start_time"][:].astype(float)
        ss_stop = ss["stop_time"][:].astype(float)
    else:
        ss_label = np.array([])
        ss_start = np.array([])
        ss_stop = np.array([])

    # --- units (hippocampus only, if a location column exists) ---
    has_location_col = "units" in f and "location" in f["units"]
    if has_location_col:
        u = f["units"]
        loc = np.array([x.decode() for x in u["location"][:]])
        hpc_idx = np.where(loc == "hpc")[0]
        sp_idx_all = u["spike_times_index"][:]
        sp_data_all = u["spike_times"]
        spike_list = []
        for i in hpc_idx:
            lo = sp_idx_all[i - 1] if i > 0 else 0
            hi = sp_idx_all[i]
            spike_list.append(np.asarray(sp_data_all[lo:hi], dtype=np.float64))
        n_units = len(spike_list)
        if n_units == 0:
            spikes_concat = np.array([], dtype=np.float64)
            spikes_bounds = np.array([0], dtype=np.int64)
        else:
            spikes_concat = np.concatenate(spike_list)
            spikes_bounds = np.concatenate([[0], np.cumsum([len(s) for s in spike_list])])
    else:
        n_units = 0
        spikes_concat = np.array([], dtype=np.float64)
        spikes_bounds = np.array([0], dtype=np.int64)

    # --- position ---
    has_position = task_group is not None
    if has_position:
        pos = f[f"processing/behavior/{task_group}/position"]
        pos_t = pos["timestamps"][:].astype(np.float64)
        pos_xy = pos["data"][:, :2].astype(np.float64)
    else:
        pos_t = np.array([])
        pos_xy = np.array([]).reshape(0, 2)

    f.close()

    np.savez_compressed(
        outfile,
        session_start_time=session_start_time,
        subject_id=subject_id,
        path=path,
        task_group=str(task_group),
        has_epochs=has_epochs,
        has_sleep_states=has_sleep_states,
        has_location_col=has_location_col,
        has_position=has_position,
        ep_labels=np.array(ep_labels),
        ep_start=ep_start,
        ep_stop=ep_stop,
        ss_label=ss_label,
        ss_start=ss_start,
        ss_stop=ss_stop,
        n_units=n_units,
        spikes_concat=spikes_concat,
        spikes_bounds=spikes_bounds,
        pos_t=pos_t,
        pos_xy=pos_xy,
    )
    usable = has_epochs and has_sleep_states and has_position and n_units >= 5
    print(f"OK {outfile} subj={subject_id} task={task_group} n_hpc_units={n_units} "
          f"n_pos={len(pos_t)} has_epochs={has_epochs} has_sleep={has_sleep_states} usable={usable}",
          flush=True)


if __name__ == "__main__":
    dandiset_id = "000213"
    path = sys.argv[1]
    outfile = sys.argv[2]
    try:
        cache_session(dandiset_id, path, outfile)
    except Exception as e:
        print(f"FAIL {path}: {e}", flush=True)
        sys.exit(1)
