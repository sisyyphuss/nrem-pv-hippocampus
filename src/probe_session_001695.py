"""
Lightweight structural + metadata probe for DANDI:001695 (Gonzalez lab)
sessions. Mirrors probe_session.py's role for 000213: run BEFORE committing
to any subject/session for the PV analysis, so that decisions ("is this a
repeated-familiar-environment session or a genuinely novel maze?") are made
from real per-file metadata, not from the dandiset-level description.

Records, per session:
  - session_description (verbatim) -- the field the task asked us to check
    per-subject/per-session rather than assuming from one sampled file.
  - whether processing/behavior/AnimalPosition/Position and .../Speed exist
    (paths may not be constant across sessions -- verify, don't assume).
  - whether intervals/SleepStates exists and its distinct `state` values.
  - units table columns and, if present, the distinct cell_type values and
    counts (Pyramidal Cell / Narrow Interneuron / Wide Interneuron per the
    task brief).
  - whether intervals/epochs exists (the task brief says this dataset likely
    lacks an explicit PRE/maze/POST epoch table -- verify per session).
"""
import sys, json
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


def probe(dandiset_id, path):
    url = get_download_url(dandiset_id, path)
    rf = remfile.File(url)
    f = h5py.File(rf, "r")
    out = {"path": path}
    out["session_start_time"] = _decode(f["session_start_time"][()]) if "session_start_time" in f else None
    out["session_description"] = _decode(f["session_description"][()]) if "session_description" in f else None
    out["subject_id"] = _decode(f["general/subject/subject_id"][()]) if "general/subject/subject_id" in f else None
    out["species"] = _decode(f["general/subject/species"][()]) if "general/subject/species" in f else None

    # --- position / speed paths: verify, do not assume constant ---
    pos_path = "processing/behavior/AnimalPosition/Position"
    spd_path = "processing/behavior/AnimalPosition/Speed"
    out["has_position_expected_path"] = pos_path in f
    out["has_speed_expected_path"] = spd_path in f
    # also record what IS under processing/behavior, in case the path differs
    if "processing/behavior" in f:
        out["behavior_module_keys"] = list(f["processing/behavior"].keys())
        sub_keys = {}
        for k in f["processing/behavior"].keys():
            try:
                sub_keys[k] = list(f[f"processing/behavior/{k}"].keys())
            except Exception:
                sub_keys[k] = None
        out["behavior_submodule_keys"] = sub_keys
    else:
        out["behavior_module_keys"] = []
        out["behavior_submodule_keys"] = {}

    if out["has_position_expected_path"]:
        pos = f[pos_path]
        out["position_n_samples"] = int(pos["data"].shape[0])
        out["position_t0"] = float(pos["timestamps"][0])
        out["position_t1"] = float(pos["timestamps"][-1])

    # --- sleep states ---
    ss_path = "intervals/SleepStates"
    out["has_sleepstates_expected_path"] = ss_path in f
    if out["has_sleepstates_expected_path"]:
        ss = f[ss_path]
        cols = list(ss.keys())
        out["sleepstates_columns"] = cols
        # NWB TimeIntervals with a custom enumerated 'state' column may be
        # ragged/vlen strings, ints, or an indexed vocabulary -- inspect the
        # raw dtype before assuming label semantics.
        if "state" in ss:
            raw = ss["state"][:]
            if raw.dtype.kind in ("S", "O"):
                vals = np.array([_decode(x) for x in raw])
            else:
                vals = raw
            uniq, counts = np.unique(vals, return_counts=True)
            out["sleepstates_state_values"] = {str(u): int(c) for u, c in zip(uniq, counts)}
        out["sleepstates_n_rows"] = int(ss["start_time"].shape[0])

    # --- epochs (explicit PRE/maze/POST table)? ---
    out["has_epochs"] = "intervals/epochs" in f
    out["intervals_keys"] = list(f["intervals"].keys()) if "intervals" in f else []

    # --- units table: cell_type column ---
    if "units" in f:
        ucols = list(f["units"].keys())
        out["units_columns"] = ucols
        out["n_units_total"] = int(f["units/id"].shape[0])
        if "cell_type" in f["units"]:
            raw = f["units/cell_type"][:]
            vals = np.array([_decode(x) for x in raw]) if raw.dtype.kind in ("S", "O") else raw
            uniq, counts = np.unique(vals, return_counts=True)
            out["cell_type_values"] = {str(u): int(c) for u, c in zip(uniq, counts)}
        if "location" in f["units"]:
            raw = f["units/location"][:]
            vals = np.array([_decode(x) for x in raw]) if raw.dtype.kind in ("S", "O") else raw
            uniq, counts = np.unique(vals, return_counts=True)
            out["unit_location_values"] = {str(u): int(c) for u, c in zip(uniq, counts)}
    else:
        out["units_columns"] = []
        out["n_units_total"] = 0

    f.close()
    return out


if __name__ == "__main__":
    try:
        result = probe("001695", sys.argv[1])
        print(json.dumps(result), flush=True)
    except Exception as e:
        print(json.dumps({"path": sys.argv[1], "error": str(e)}), flush=True)
