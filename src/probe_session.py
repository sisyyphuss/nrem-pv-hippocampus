"""
Lightweight structural probe (metadata only, no big array reads) of a
DANDI:000213 session: which behavioral task module it uses, how many
hippocampal units it has, and whether before/during/after-trial epochs
and sleep_states exist. Used for QC/session-selection BEFORE committing
to a subject for the main analysis.
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
    asset_id = d["results"][0]["asset_id"]
    r2 = requests.get(
        f"https://api.dandiarchive.org/api/dandisets/{dandiset_id}/versions/{version}/assets/{asset_id}/",
        timeout=60,
    )
    r2.raise_for_status()
    urls = r2.json()["contentUrl"]
    return urls[1] if len(urls) > 1 else urls[0]


def probe(dandiset_id, path):
    url = get_download_url(dandiset_id, path)
    rf = remfile.File(url)
    f = h5py.File(rf, "r")
    out = {"path": path}
    out["session_start_time"] = f["session_start_time"][()].decode()
    out["subject_id"] = f["general/subject/subject_id"][()].decode()

    beh_keys = list(f["processing/behavior"].keys()) if "processing/behavior" in f else []
    task_key = [k for k in beh_keys if k != "allocentric_frame_tracking"]
    out["task_group"] = task_key[0] if task_key else None
    out["has_position"] = bool(task_key and "position" in f[f"processing/behavior/{task_key[0]}"])

    out["has_sleep_states"] = "processing/ecephys/sleep_states" in f
    out["has_epochs"] = "intervals/epochs" in f
    if out["has_epochs"]:
        ep = f["intervals/epochs"]
        tags_flat = ep["tags"][:]
        tags_idx = ep["tags_index"][:]
        labels = []
        prev = 0
        for i in range(len(tags_idx)):
            tg = tags_flat[prev:tags_idx[i]]
            prev = tags_idx[i]
            labels.append(tg[0].decode() if len(tg) else "")
        out["epoch_labels"] = labels
    else:
        out["epoch_labels"] = []

    if "units" in f and "location" in f["units"]:
        loc = np.array([x.decode() for x in f["units/location"][:]])
        vals, counts = np.unique(loc, return_counts=True)
        out["unit_locations"] = dict(zip(vals.tolist(), [int(c) for c in counts]))
        out["n_hpc_units"] = int(out["unit_locations"].get("hpc", 0))
    else:
        out["unit_locations"] = {}
        out["n_hpc_units"] = 0 if "units" not in f else int(f["units/id"].shape[0])
        out["n_units_no_location_col"] = int(f["units/id"].shape[0]) if "units" in f else 0

    f.close()
    print(json.dumps(out), flush=True)


if __name__ == "__main__":
    try:
        probe("000213", sys.argv[1])
    except Exception as e:
        print(json.dumps({"path": sys.argv[1], "error": str(e)}), flush=True)
