"""
Explicit, chronologically ordered session paths for DANDI:001695 (Gonzalez
lab), retained for Analysis B / C, plus the exclusion reasoning that was
actually verified against per-file metadata (not the dandiset-level
description alone -- see notebook Analysis B, "Novel maze" investigation).

--------------------------------------------------------------------------
WHAT WAS PROBED (all 22 assets in the dandiset, via probe_session_001695.py,
results in probe_results_001695.jsonl):

Subjects and session TYPES found:
  - M01, M02, M03, M05: each has a mix of
      (a) "Head Fix Session with Odor Stimulation" (`*_ecephys.nwb`, no
          `processing/behavior` module at all -- no position, no sleep
          scoring). EXCLUDED: not usable for this analysis at all (no
          spatial code to correlate, no sleep to quantify).
      (b) "Novel Maze Session with Pre and Post sleep" (`*_behavior+ecephys.nwb`)
          -- has position (`processing/behavior/AnimalPosition/Position`),
          a separate `processing/behavior/Speed` series, and
          `intervals/SleepStates`.
  - NN-syn-M02, NN-syn-M05: "Head-fix SiNAPS Recording" only, no
    `processing/behavior` module at all. EXCLUDED entirely (no position, no
    sleep scoring, not even usable for the 000061-style descriptive-only
    role since that role is already filled by DANDI:000061 per the task
    scope).

--------------------------------------------------------------------------
THE "NOVEL MAZE" QUESTION (the task's explicit ask: verify precisely,
don't assume from one sampled file):

EVERY SINGLE `behavior+ecephys` session, for EVERY subject (M01, M02, M03,
M05), has the IDENTICAL `session_description` string:
  "Novel Maze Session with Pre and Post sleep"
This text is therefore UNINFORMATIVE for distinguishing "same environment
reused" from "genuinely new maze each day" -- it does not vary at all
across days, so it cannot be used as the evidence source the task asked
for ("document your exact evidence and reasoning per subject").

The dandiset-level overview (fetched via the DANDI API `description` field)
states: "Neuropixels sessions include PRE-rest, Maze navigation (Novel or
2-maze configurations), and POST-rest homecage recordings." -- i.e. the
dandiset explicitly acknowledges BOTH a novel-maze design AND a distinct
"2-maze" (repeated, alternating) design, but does not say which files are
which.

RESOLUTION ADOPTED HERE: since the text metadata cannot disambiguate, the
actual `processing/behavior/AnimalPosition/Position` (x, y) data range
(footprint) for every maze session was pulled directly and compared within
each subject. The apparatus track length (x-span) is consistently
~97-117 raw position units on every session of every subject (i.e. one
fixed-length linear/linearized track is used throughout for M01/M02/M03/M05
-- this is NOT "a genuinely different maze/apparatus every day"). But the
TRACK WIDTH (y-span, and equivalently the aspect ratio y_span/x_span) is
clearly BIMODAL within-subject:

  subject  date        x_span  y_span  aspect(y/x)  cluster
  M01      2024-03-12   111.1    7.7    0.069        narrow
  M01      2024-03-13   108.7    5.4    0.050        narrow
  M01      2024-03-14   109.5    5.7    0.052        narrow
  M01      2024-03-18   108.0   10.3    0.095        WIDE   <-- config change
  M02      2024-03-12   107.7    6.0    0.056        narrow
  M02      2024-03-13   108.7    5.4    0.050        narrow
  M02      2024-03-14   106.7    5.2    0.049        narrow
  M02      2024-03-18   110.6   10.5    0.095        WIDE   <-- config change
  M03      2024-06-21   111.7    5.8    0.052        narrow
  M03      2024-06-22    97.5   14.2    0.146        WIDE   <-- alternates
  M03      2024-06-23   111.5    8.4    0.075        narrow    every
  M03      2024-06-24   108.0   13.1    0.121        WIDE   <-- other day
  M05      2024-07-29     1.1    0.1    (excluded -- anomalous position
                                          scale, ~100x smaller raw span than
                                          every other session in the whole
                                          dandiset; SpatialSeries `unit`/
                                          `conversion` metadata attributes
                                          are IDENTICAL to the other M05
                                          files (unit="centimeters",
                                          conversion=0.01), so this is not a
                                          units-label artifact -- it looks
                                          like a real tracking/calibration
                                          failure for this specific file.
                                          Excluded as a QC failure, not an
                                          environment-identity judgement.)
  M05      2024-07-30   100.0   12.8    0.128        WIDE
  M05      2024-07-31   112.4    7.5    0.067        narrow

A clean gap exists in the aspect-ratio distribution between 0.075 (M03
6/23, narrow) and 0.095 (M01/M02 3/18, WIDE) across all 13 usable sessions,
so `aspect_ratio < 0.085` vs `>= 0.085` is used as an explicit, data-driven
(not literature-borrowed) same-environment classifier. This gives, per
subject, exactly the "2-maze configuration" pattern the dandiset overview
warns about:
  - M01: narrow narrow narrow WIDE  -> 3/18 is a different configuration.
  - M02: narrow narrow narrow WIDE  -> 3/18 is a different configuration
         (both M01 and M02 change on the SAME calendar day -- consistent
         with a cohort-wide apparatus change on 3/18, not an animal-
         specific effect).
  - M03: narrow WIDE narrow WIDE    -> alternates every session; this looks
         exactly like the "2-maze configuration" the dandiset overview
         describes, i.e. genuinely different physical arrangements on
         interleaved days for this subject.
  - M05: (7/29 excluded) WIDE narrow -> the only two usable sessions are in
         DIFFERENT configurations, so M05 contributes ZERO same-environment
         pairs.

CONSEQUENCE FOR PAIR-BUILDING (implemented in build_pairs.py,
`build_pair_table_001695()`): sessions are grouped, per subject, into
same-`aspect_ratio` clusters (preserving chronological order within each
cluster), and only CONSECUTIVE sessions WITHIN the same cluster are paired
-- exactly analogous to how `build_pair_table()` already restricts
DANDI:000213 pairs to the subject's dominant task label, and allows the
day-gap to exceed 1 when a differently-labelled day is skipped. This is a
real methodological weak point relative to DANDI:000213 (which has an
explicit before/during/after-trials epoch table per day, not a
geometry-inferred proxy for "same environment") and is disclosed as such in
the notebook, not hidden.

--------------------------------------------------------------------------
SUBJECT/SESSION INCLUSION TABLE (mirrors the 000213 QC table style):

  subject | n_ecephys-only (head-fix, excluded) | n_maze_sessions | n_usable_maze | reason for any exclusion
  M01     | 1                                    | 4               | 4             | none -- all 4 usable, but 3/18 forms its own 1-session cluster (no pair partner within-cluster beyond the narrow trio)
  M02     | 2                                    | 4               | 4             | same as M01
  M03     | 0                                    | 4               | 4             | alternating clusters -- every session usable, but only same-cluster consecutive pairs qualify
  M05     | 1                                    | 3               | 2             | 2024-07-29 excluded: anomalous position-tracking scale (QC failure, see above)
  NN-syn-M02 | 1 (all sessions)                  | 0               | 0             | no `processing/behavior` module at all -- excluded entirely
  NN-syn-M05 | 2 (all sessions)                  | 0               | 0             | same as NN-syn-M02
"""

# Only behavior+ecephys sessions (position + sleep present) are listed here.
# Head-fix odor-stimulation and SiNAPS-only sessions are intentionally
# omitted (see module docstring).
SESSIONS_001695 = {
    "M01": [
        "sub-M01/sub-M01_ses-20240312T100000_behavior+ecephys.nwb",
        "sub-M01/sub-M01_ses-20240313T100000_behavior+ecephys.nwb",
        "sub-M01/sub-M01_ses-20240314T100000_behavior+ecephys.nwb",
        "sub-M01/sub-M01_ses-20240318T100000_behavior+ecephys.nwb",
    ],
    "M02": [
        "sub-M02/sub-M02_ses-20240312T100000_behavior+ecephys.nwb",
        "sub-M02/sub-M02_ses-20240313T100000_behavior+ecephys.nwb",
        "sub-M02/sub-M02_ses-20240314T100000_behavior+ecephys.nwb",
        "sub-M02/sub-M02_ses-20240318T100000_behavior+ecephys.nwb",
    ],
    "M03": [
        "sub-M03/sub-M03_ses-20240621T100000_behavior+ecephys.nwb",
        "sub-M03/sub-M03_ses-20240622T100000_behavior+ecephys.nwb",
        "sub-M03/sub-M03_ses-20240623T100000_behavior+ecephys.nwb",
        "sub-M03/sub-M03_ses-20240624T100000_behavior+ecephys.nwb",
    ],
    "M05": [
        # 2024-07-29 intentionally omitted -- anomalous position-tracking
        # scale (QC failure), see module docstring.
        "sub-M05/sub-M05_ses-20240730T100000_behavior+ecephys.nwb",
        "sub-M05/sub-M05_ses-20240731T100000_behavior+ecephys.nwb",
    ],
}

# Same-environment aspect-ratio classifier threshold, derived from the
# observed gap in the data (0.075 vs 0.095) -- see module docstring.
ASPECT_RATIO_THRESHOLD = 0.085
