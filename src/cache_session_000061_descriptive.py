"""
DESCRIPTIVE-ONLY cache of DANDI:000061 (Girardeau lab) sessions.

Per explicit user instruction, DANDI:000061 is NOT used for the PV-correlation
analysis: it ships raw video instead of extracted x,y position, so place-field
rate maps (and therefore any PV correlation) cannot be computed without a
separate, substantial video-tracking step that is out of scope here. This
dataset is used ONLY for a descriptive/QC section: recording-day counts per
rat, and a sleep-architecture summary (NREM/REM/Wake durations per day) built
from `processing/behavior/states`.

IMPORTANT UNIT-MISMATCH FINDING (investigated before writing any downstream
number from this dataset -- see notebook Descriptive-only section for the
full writeup):

  `intervals/epochs` (start_time/stop_time) and `processing/behavior/states`
  (start_time/stop_time) do NOT use the same time units in this dataset.
  Verified directly (not assumed) on sub-Rat08/.../20130708:
    - intervals/epochs stop_time max  = 24,759,200  (raw)
    - processing/behavior/states stop_time max = 24,754.0  (raw)
    - 24,759,200 / 1000 = 24,759.2 -- matches the states-table max (24,754.0)
      to within 5 seconds (a real, small, boundary rounding difference, not
      noise at the 1000x scale).
    - intervals/epochs tags are 'sleep'/'run' (a coarse 2-way behavioral
      split), each with only ~2-4 rows per session -- consistent with a
      block-epoch table given in MILLISECONDS.
    - processing/behavior/states has Awake/Non-REM(/REM) rows on the order of
      tens-to-hundreds of seconds each, with start times as small as 5.0 --
      consistent with a fine-grained state table given in SECONDS.
  Resolution adopted here: `intervals/epochs` start_time/stop_time are
  divided by 1000 before being combined with anything from `states`. The
  `states` table is treated as already being in seconds (it is used as-is).
  This cache stores BOTH tables in their ORIGINAL raw units (no conversion
  applied at cache time) plus an explicit `epochs_time_unit_correction`
  scalar (1000.0) so the conversion is applied explicitly, visibly, and only
  once, in the notebook -- never silently baked into a cached file.

Usage: python3 cache_session_000061_descriptive.py <path_in_dandiset> <outfile.npz>
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


def cache_session(dandiset_id, path, outfile):
    url = get_download_url(dandiset_id, path)
    rf = remfile.File(url)
    f = h5py.File(rf, "r")

    session_start_time = _decode(f["session_start_time"][()]) if "session_start_time" in f else ""
    subject_id = _decode(f["general/subject/subject_id"][()]) if "general/subject/subject_id" in f else ""

    # --- processing/behavior/states: Awake / Non-REM / REM epochs, SECONDS ---
    has_states = "processing/behavior/states" in f
    if has_states:
        st = f["processing/behavior/states"]
        st_label = np.array([_decode(x) for x in st["label"][:]])
        st_start = st["start_time"][:].astype(float)
        st_stop = st["stop_time"][:].astype(float)
    else:
        st_label, st_start, st_stop = np.array([]), np.array([]), np.array([])

    # --- intervals/epochs: coarse sleep/run blocks, MILLISECONDS (raw, unconverted) ---
    has_epochs = "intervals/epochs" in f
    ep_labels, ep_start, ep_stop = [], np.array([]), np.array([])
    if has_epochs:
        ep = f["intervals/epochs"]
        tags_flat = ep["tags"][:]
        tags_idx = ep["tags_index"][:]
        ep_start = ep["start_time"][:].astype(float)
        ep_stop = ep["stop_time"][:].astype(float)
        prev = 0
        for i in range(len(ep_start)):
            tg = tags_flat[prev:tags_idx[i]]
            prev = tags_idx[i]
            ep_labels.append(_decode(tg[0]) if len(tg) else "")

    # --- units: hpc unit count only, for QC context (not used for any PV computation) ---
    n_hpc_units = 0
    if "units" in f and "location" in f["units"]:
        loc = np.array([_decode(x) for x in f["units/location"][:]])
        n_hpc_units = int((loc == "hpc").sum())

    f.close()

    np.savez_compressed(
        outfile,
        session_start_time=session_start_time,
        subject_id=subject_id,
        path=path,
        has_states=has_states,
        st_label=st_label, st_start=st_start, st_stop=st_stop,
        has_epochs=has_epochs,
        ep_labels=np.array(ep_labels), ep_start=ep_start, ep_stop=ep_stop,
        epochs_time_unit_correction=1000.0,  # divide ep_start/ep_stop by this to get seconds
        n_hpc_units=n_hpc_units,
    )
    print(f"OK {outfile} subj={subject_id} n_states_rows={len(st_label)} "
          f"n_epoch_rows={len(ep_labels)} n_hpc_units={n_hpc_units}", flush=True)


if __name__ == "__main__":
    dandiset_id = "000061"
    path = sys.argv[1]
    outfile = sys.argv[2]
    try:
        cache_session(dandiset_id, path, outfile)
    except Exception as e:
        print(f"FAIL {path}: {e}", flush=True)
        sys.exit(1)
