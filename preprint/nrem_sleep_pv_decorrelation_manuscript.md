Testing whether inter-session NREM sleep predicts hippocampal population-vector decorrelation: an open, multi-dataset, multi-species replication attempt
=====================================================================================================================================================

**Saad Yusuf**  
*Department of Computer Science, Jamia Millia Islamia, New Delhi, India*  
Correspondence: `syed.saad7527@gmail.com` | ORCID: `0009-0000-2333-0078`

*Preprint — not yet peer reviewed.*

---

## Abstract

Two recent studies independently reported that non-rapid-eye-movement (NREM) sleep destabilizes, or "decorrelates," the hippocampal population code for space: Valdivia et al. (2024) found that a population-vector-based spatial stability index correlated negatively with intervening NREM duration in rats performing an object-place recognition task (r = −0.55, p = 0.028), and Bollmann et al. (2025) reported that NREM sleep accelerates hippocampal assembly-pattern "drift" while REM sleep opposes it. Both studies used purpose-collected recordings with a single, short (90-minute to ~20-hour) intervening rest period. We asked whether the same relationship is detectable in openly archived, chronic, multi-day hippocampal recordings, using data not collected for this question. We built an analysis pipeline around two public DANDI Archive dandisets — DANDI:000213 (rat dorsal CA1 and lateral septum, Tingley et al., Buzsáki Lab, 3 subjects retained after quality control, 20 same-task session-pairs) and DANDI:001695 (mouse CA1, Gonzalez et al., 3 subjects retained, 6 same-configuration session-pairs) — computing an ensemble-tuning-curve population-vector (PV) correlation between same-environment session pairs and regressing it against the NREM duration recorded between those sessions. Across four pre-specified analyses (each dataset alone, the two datasets combined with dataset/species as a covariate, and a pyramidal-versus-interneuron comparison within DANDI:001695), we found a consistent negative point estimate (combined Spearman r = −0.246, p = 0.225, n = 26 pairs, 6 subjects) but no result that survived a battery of robustness checks: within-subject permutation testing, an independent identity-free similarity metric (the RV-coefficient), leave-one-subject/dataset-out re-fitting, and restriction to consecutive-day pairs. We report this as an inconclusive, underpowered replication attempt rather than a confirmatory or disconfirmatory finding, and we make the full pipeline — including explicit documentation of every methodological compromise required by data not originally collected for this purpose — openly available as a reusable resource for adequately powered follow-up work.

**Keywords:** hippocampus, CA1, sleep, NREM, representational drift, population vector correlation, open data, DANDI, replication

---

## 1. Introduction

The hippocampal population code for space is not static: place fields drift, reorganize, and partially "remap" across repeated exposures to the same environment, even in the absence of any experimentally imposed change (Ziv et al., 2013). A growing body of work implicates sleep — and NREM sleep specifically — as an active driver of this reorganization rather than a passive consolidation window that simply preserves the code unchanged. Most directly relevant to the present study, Valdivia et al. (2024) recorded dorsal CA1 in rats performing an object-place recognition task bracketed by a single ~90-minute sleep opportunity and found that a population-vector-based stability index was negatively correlated with the amount of NREM obtained (r = −0.55, p = 0.028); notably, this decorrelation was present only when the environment configuration changed between encoding and retrieval, and was absent in a stable open-field control condition. Bollmann et al. (2025) recorded CA1 across ~17–20 hours of post-learning rest and found that reactivation of previously learned firing patterns gradually gave way to patterns resembling future, not-yet-experienced sessions, and that this drift was accelerated specifically during NREM epochs and opposed during REM epochs.

Both studies are compelling but share a common design: a single learning episode bracketed by one, relatively short, intervening rest/sleep period, recorded in one species (rat), in one lab. This leaves open a related but distinct question: in animals that repeatedly re-visit the same familiar environment across many days — the classic design used to study place-field stability and representational drift more broadly (Ziv et al., 2013) — does the amount of NREM sleep recorded *between* any two such visits predict how much the population code has reorganized by the next visit? This is a within-animal, repeated-measures version of the same hypothesis, and it is, in principle, testable using chronic multi-day recordings that already exist in public archives, without new data collection.

We used the DANDI Archive (Neurodata Without Borders, NWB, format) to identify and rigorously vet candidate open datasets against a specific set of requirements: chronic hippocampal recordings from the same animals across multiple days in a stable, repeated environment; simultaneously recorded, spike-sorted units sufficient to build a population code; continuous position tracking; scored NREM/REM/wake sleep states; and — critically — verified, real recording timestamps rather than inferred or placeholder values, since the entire predictor variable (time spent in NREM *between* two specific recordings) depends on trusting those timestamps. We deliberately did not restrict this search to datasets already associated with a sleep/drift publication, in order to test the hypothesis on data collected for other purposes.

This paper reports what we found: a real, working, cross-validated analysis pipeline, applied to two independent open datasets from two labs and two rodent species, and an honest null/inconclusive result. We frame this explicitly as (i) a replication attempt of Valdivia et al. (2024) and Bollmann et al. (2025) in a different experimental regime (chronic, many-repetition, multi-day exposure rather than a single learning episode), (ii) a methodological resource — an open, documented, reusable pipeline for this class of question — and (iii) pilot data characterizing the effect size and variance that a properly powered follow-up study would need to plan around.

---

## 2. Methods

### 2.1 Datasets and inclusion criteria

We required candidate datasets to provide, verified at the level of the raw data files rather than the dataset description: (a) chronic extracellular hippocampal recordings from the same subject across ≥5 recording days; (b) a spike-sorted units table; (c) continuous position tracking; (d) either pre-scored NREM/REM/wake sleep states or sufficient LFP/EMG to score them; (e) real, verifiable recording timestamps (not placeholder or inferred dates); and (f) a behavioral paradigm involving repeated exploration of a nominally constant, familiar environment. Two DANDI Archive dandisets satisfied all criteria and are used for the primary analyses below:

- **DANDI:000213** (Tingley et al., Buzsáki Laboratory), rat, silicon-probe recordings from dorsal CA1 and lateral septum, with an explicit `before trials` / `during trials` / `after trials` epoch structure recorded every day (i.e., sleep/rest recorded on both sides of each day's behavioral session) and a dedicated sleep-state interval table (Awake / Non-REM / REM).
- **DANDI:001695** (Gonzalez et al.), mouse, high-density Neuropixels/SiNAPS recordings from CA1 (and, in other sessions not used here, DG/CA2/CA3/retrosplenial cortex), with a `cell_type`-annotated units table (Pyramidal Cell / Narrow Interneuron / Wide Interneuron) and a sleep-state interval table (NREM / REM / WAKE / Ripple).

A third dandiset, **DANDI:000061** (Girardeau et al., Buzsáki Laboratory; hippocampus and basolateral amygdala, 4 rats, 7–13 consecutive daily sessions each, with an explicit sleep–run–sleep–run epoch structure on some days), met every criterion above except (c): position was recorded only as raw video, with no extracted coordinate trace. Because pose extraction from video was outside the scope of this analysis, DANDI:000061 was retained only for a descriptive summary of its sleep architecture (Section 3.0) and did not contribute any session-pairs to the population-vector analyses.

### 2.2 Data quality control and session/subject selection

Every session file in both primary dandisets was probed programmatically (via the DANDI API and direct, byte-range NWB streaming with `remfile`/`h5py`) before any subject or session was included, rather than relying on dataset-level descriptions. For DANDI:000213, this excluded two of five subjects: one had no sleep-state scoring in any of its 17 sessions, and one had a task apparatus that rotated across up to four different physical configurations with degraded hippocampal unit yield in its later sessions. The three retained subjects (11, 16, and 9 sessions respectively) were further restricted, per subject, to sessions sharing that subject's dominant task label.

For DANDI:001695, the session-level metadata field intended to describe the environment (`session_description`) was found to be an identical, uninformative string ("Novel Maze Session with Pre and Post sleep") for every session of every subject, and so could not be used to determine which sessions shared a physical environment. We instead used the raw position-tracking footprint (the spatial extent of recorded x,y coordinates) as direct evidence: track length was consistently ~97–117 position units across all sessions of all subjects (i.e., a single fixed-length track geometry throughout), but track width was clearly bimodal within each subject, revealing the two-configuration ("2-maze") design referenced in the dandiset's own overview text without being tied to specific session identifiers. A data-driven aspect-ratio threshold, chosen from a clear gap in the observed distribution, was used to cluster each subject's sessions into same-configuration groups; only consecutive sessions within the same cluster were treated as a same-environment pair. One session (subject M05) was excluded for an anomalous, ~100-fold smaller position-tracking scale consistent with a calibration failure rather than a genuine change in tracking units (confirmed by identical `SpatialSeries` unit/conversion metadata to that subject's other sessions). This yielded 6 usable same-environment, consecutive-day pairs across 3 subjects (M01, M02, M03); the fourth retained-eligible subject (M05) contributed zero pairs because its only two artifact-free sessions were in different configurations.

### 2.3 Rate maps and population-vector correlation

For each session, position was binned into a shared two-dimensional grid whose edges were fixed from the pooled position range of the specific pair being compared (so that bins align across the two days being correlated). For every hippocampal unit, an occupancy-normalized firing-rate map was constructed from spike counts divided by occupancy time within the behavioral epoch (DANDI:000213: the explicit `during trials` window; DANDI:001695, lacking an explicit epoch table, the window in which position tracking was actively recorded, treated as a proxy for the maze bout). Spatial bins with insufficient sampling (<0.15 s of occupancy) in either session of a pair were excluded from that pair's comparison.

Because neither dataset supports tracking individual cell identity across recording days — DANDI:000213's silicon probe is advanced in depth across sessions, and no equivalent identity information is available for DANDI:001695 — a literal, cell-identity-matched population vector correlation (in the sense used by, e.g., Leutgeb et al., 2005) could not be constructed. We instead computed, for each session, an **ensemble tuning curve**: each included unit's rate map was z-scored across valid spatial bins (removing across-unit differences in overall firing rate, which can also differ systematically across days due to electrode drift) and then averaged across all included units, yielding a single one-dimensional vector per session (length = number of spatial bins valid in both sessions of the pair). Units with mean firing rate below 0.1 Hz in the epoch were excluded as too sparse to inform a rate map; for DANDI:001695, units were additionally pre-filtered by their experimenter-assigned `cell_type` label (Pyramidal Cell for the primary analysis; Narrow/Wide Interneuron for the interneuron comparison in Section 3.4) rather than by a firing-rate heuristic, since DANDI:000213 lacks any equivalent cell-type annotation. The Pearson correlation between the two sessions' ensemble tuning curves is the primary outcome variable, referred to throughout as the PV correlation.

As an independent, methodologically distinct robustness check that does not require collapsing cells into a single average and does not require matched cell counts between sessions, we additionally computed the **RV-coefficient** (Robert & Escoufier, 1976) — a multivariate generalization of squared correlation — between each pair's full [bins × cells] rate-map matrices.

### 2.4 Quantifying NREM sleep between sessions

For a pair of sessions (day N, day N+1 in the subject's retained, same-environment sequence), NREM duration "between" the two recordings was calculated from each dataset's sleep-state interval table as: NREM recorded during the `after trials` (DANDI:000213) or post-tracking (DANDI:001695) window of day N, plus the full duration of any NREM recorded during any strictly intervening recorded day (present for some but not all pairs), plus NREM recorded during the `before trials` / pre-tracking window of day N+1. This is necessarily a lower-bound proxy: any sleep occurring in the animal's home cage outside the NWB-recorded windows is invisible to the archived data and cannot be recovered. REM and wake durations were computed identically for use as covariates.

### 2.5 Statistical analysis

The primary test in every analysis was the Spearman rank correlation between between-session NREM duration and PV correlation, pooled across all usable pairs (Pearson correlation is reported alongside as a secondary check). Because multiple pairs are contributed by the same subject and are therefore not independent observations, we additionally ran: (i) a within-subject permutation test (NREM values shuffled within each subject's own set of pairs, 5,000 iterations, preserving each subject's marginal distributions of both variables while destroying their pairing); (ii) a linear mixed-effects model with subject as a random intercept and wake duration, REM duration, and calendar-day gap as fixed-effect covariates, to test whether any NREM effect was specific rather than a proxy for elapsed time or total sleep; (iii) leave-one-subject-out and, for the combined analysis, leave-one-dataset-out re-fitting; (iv) the RV-coefficient described above, correlated against NREM duration independently of the primary ensemble-tuning-curve metric; and (v) restriction to the subset of pairs separated by exactly one calendar day, to reduce the influence of any unmonitored intervening sleep. For the pyramidal-versus-interneuron comparison (Section 3.4), paired differences in PV correlation were additionally summarized by sign count and, given the sample size, a Wilcoxon signed-rank test reported descriptively rather than as a confirmatory significance test. No correction for multiple comparisons across this battery is applied; the battery is reported in full, and results should be read as a robustness profile rather than as five independent chances to find significance.

All analyses were implemented in Python (numpy, pandas, scipy, statsmodels) and are available in full, including every intermediate quality-control decision, in the accompanying Jupyter notebook (Section 6).

---

## 3. Results

### 3.0 DANDI:000061 — descriptive sleep architecture (not used in the PV analysis)

Across the four DANDI:000061 subjects retained after quality control, one subject (Rat09) had usable sleep-state scoring in only 1 of its 10 sessions; the remaining three subjects had complete or near-complete sleep-state coverage across 7–13 consecutive daily sessions each, with a clear NREM/wake/REM architecture recorded on both sides of most daily behavioral sessions, and, on some days, an explicit sleep epoch recorded *between* two behavioral run epochs within the same day. This dataset is structurally the closest match in the wider DANDI archive to the literal framing of "NREM between two recording sessions," and is flagged as a priority target for a future extension of this pipeline once position tracking is available (Section 5).

![Figure 1](../figures/09_000061_sleep_architecture.png)
**Figure 1. DANDI:000061 sleep architecture (descriptive only).** Per-subject, per-day summary of recorded Awake / Non-REM / REM duration across the four retained DANDI:000061 subjects, illustrating that sleep-state coverage is near-complete for three of the four subjects and that a sleep-run-sleep-run epoch structure is present on a subset of days. This dataset did not contribute session-pairs to Analyses A–D because it lacks extracted position tracking (Section 2.1).

### 3.1 Analysis A — DANDI:000213 alone

Across 20 usable same-task session-pairs from 3 rats (DT2, DT5, DT7), the pooled correlation between NREM duration and PV correlation was small and not significant (Spearman r = −0.168, p = 0.478; Pearson r = −0.152, p = 0.522). The within-subject permutation test gave an empirical p-value of 0.468. A mixed-effects model including wake duration, REM duration, and day-gap as covariates did not identify a reliable NREM-specific effect (coefficient ≈ 0, p = 0.82) and showed convergence difficulties consistent with too few observations for the number of fixed effects. Leave-one-subject-out re-fitting produced correlations ranging from −0.31 to +0.08 depending on which single animal was held out. The independent RV-coefficient check found essentially no relationship with NREM duration (r = 0.018, p = 0.94) and did not agree in sign with the primary metric. Restricting to the 8 pairs separated by exactly one calendar day reversed the sign of the point estimate (r = +0.31, p = 0.46).

![Figure 2](../figures/01_example_hypnograms.png)
**Figure 2. Example hypnograms, DANDI:000213.** Awake / Non-REM / REM state timelines for one example same-task session pair (day N, day N+1) from subject DT2, illustrating the sleep-state interval data underlying the NREM-between-sessions calculation (Section 2.4).

![Figure 3](../figures/02_example_rate_maps.png)
**Figure 3. Example CA1 rate maps, DANDI:000213.** Occupancy-normalized firing-rate maps for a sample of hippocampal units from one example session, illustrating the spatial tuning that feeds into the ensemble tuning curve (Section 2.3).

![Figure 4](../figures/03_example_ensemble_tuning.png)
**Figure 4. Example ensemble tuning curves and pairwise PV correlation.** The two sessions' ensemble tuning curves (per-unit z-scored rate maps averaged across all included units) for the example pair in Figures 2–3, with the resulting Pearson PV correlation annotated.

![Figure 5](../figures/04_main_scatter.png)
**Figure 5. Analysis A primary result.** NREM duration recorded between sessions (x-axis) versus PV correlation (y-axis) for all 20 usable DANDI:000213 pairs, colored by subject (DT2 n=9, DT5 n=8, DT7 n=3 pairs). Spearman r = −0.168, p = 0.478.

![Figure 6](../figures/05_permutation_null.png)
**Figure 6. Within-subject permutation null distribution, Analysis A.** Null distribution of the Spearman correlation under 5,000 within-subject shuffles of NREM duration, with the observed statistic marked; empirical p = 0.468.

![Figure 7](../figures/06_leave_one_out.png)
**Figure 7. Leave-one-subject-out stability, Analysis A.** Pooled Spearman correlation recomputed after excluding each subject in turn, showing the point estimate ranges from −0.31 to +0.08 depending on which single animal is held out.

![Figure 8](../figures/07_method_robustness.png)
**Figure 8. Independent RV-coefficient cross-check, Analysis A.** Comparison of the primary ensemble-tuning-curve PV correlation against the identity-free RV-coefficient metric for the same pairs, and the RV-coefficient's own correlation with NREM duration (r = 0.018, p = 0.94).

### 3.2 Analysis B — DANDI:001695 alone

Across the 6 usable same-configuration, consecutive-day pairs from 3 mice (M01, M02, M03), the pooled Spearman correlation was numerically large (r = −0.829, p = 0.042). We flag this result explicitly as unreliable in isolation: with n = 6 pairs from 3 subjects, a single pooled correlation of this kind is exactly the pattern expected to arise by chance in a small sample and to fail to replicate; it should not be read as independent confirmation of the hypothesis.

![Figure 9](../figures/08_001695_example_hypnogram.png)
**Figure 9. Example hypnogram, DANDI:001695.** Sleep-state timeline (NREM / REM / WAKE) for one example DANDI:001695 session, illustrating the `intervals/SleepStates` data used identically to Figure 2 for this dataset.

![Figure 10](../figures/08b_analysisB_scatter.png)
**Figure 10. Analysis B primary result.** NREM duration between sessions versus pyramidal-ensemble PV correlation for all 6 usable DANDI:001695 pairs (M01, M02, M03; 2 pairs each). Spearman r = −0.829, p = 0.042 — reported with the explicit caveat above that n = 6 is too small to trust this point estimate in isolation.

### 3.3 Analysis C — Combined multi-dataset analysis

Pooling both datasets (26 pairs, 6 subjects, 2 species) gave a combined Spearman correlation of r = −0.246 (p = 0.225; Pearson r = −0.173, p = 0.399) — numerically larger in magnitude and nominally lower in p-value than DANDI:000213 alone, but still not statistically significant. The within-subject, dataset-block-preserving permutation test gave an empirical p = 0.224. Leave-one-subject-out correlations ranged from −0.42 to −0.11; leave-one-dataset-out re-fitting exactly reproduced each dataset's own single-dataset result (removing DANDI:000213 leaves DANDI:001695's r = −0.83; removing DANDI:001695 leaves DANDI:000213's r = −0.17), confirming that the combined estimate is a genuine pool of the two rather than an artifact of the merge. The independent RV-coefficient check again found essentially no relationship (r = −0.006, p = 0.975) and did not agree in sign with the primary metric. Restricting to the 12 combined pairs separated by exactly one day gave a near-zero correlation (r = −0.056, p = 0.863). Because dataset identity is confounded 1:1 with species (every DANDI:000213 subject is rat, every DANDI:001695 subject is mouse), we treat any dataset/species difference as an explicit covariate and interaction term rather than assuming species-invariance; the available sample size does not permit a reliable estimate of that interaction.

![Figure 11](../figures/10_correlation_by_dataset.png)
**Figure 11. Per-dataset versus combined pooled correlation.** Spearman r between NREM duration and PV correlation computed separately for DANDI:000213 alone, DANDI:001695 alone, and the two combined, showing the combined estimate sits between the two single-dataset estimates.

![Figure 12](../figures/11_combined_scatter.png)
**Figure 12. Analysis C primary result.** NREM duration between sessions versus PV correlation for all 26 combined pairs, colored by subject and dataset (circles: DANDI:000213/rat; triangles: DANDI:001695/mouse). Spearman r = −0.246, p = 0.225.

![Figure 13](../figures/12_combined_permutation.png)
**Figure 13. Within-subject, dataset-block-preserving permutation null, Analysis C.** Null distribution under 5,000 shuffles preserving subject and dataset structure; empirical p = 0.224.

![Figure 14](../figures/13_combined_loso_subject.png)
**Figure 14. Leave-one-subject-out stability, Analysis C.** Pooled Spearman correlation recomputed after excluding each of the 6 subjects in turn; range −0.42 to −0.11.

![Figure 15](../figures/14_combined_lodo.png)
**Figure 15. Leave-one-dataset-out stability, Analysis C.** Pooled Spearman correlation recomputed after excluding each dataset in turn; removing DANDI:000213 leaves DANDI:001695's own r = −0.83, and removing DANDI:001695 leaves DANDI:000213's own r = −0.17, confirming the combined estimate is a genuine pool of the two single-dataset results.

![Figure 16](../figures/15_combined_method_robustness.png)
**Figure 16. Independent RV-coefficient cross-check, Analysis C.** As Figure 8, computed on the combined 26-pair dataset (r = −0.006, p = 0.975).

### 3.4 Analysis D — Interneuron versus pyramidal ensembles (DANDI:001695 only)

Because only DANDI:001695 provides a waveform-based cell-type classification, we repeated the ensemble-tuning-curve PV correlation for the same 6 session-pairs using interneuron units (Narrow + Wide Interneuron) in place of pyramidal units. Raw interneuron counts per session ranged from 26 to 224, well above the minimum required by the ensemble-tuning-curve method, and all 6 pairs remained usable. The interneuron-ensemble correlation with NREM duration was weaker than the pyramidal result on the identical pairs (Spearman r = −0.429, p = 0.397, versus r = −0.829, p = 0.042 for pyramidal units). Pyramidal PV correlation exceeded interneuron PV correlation in 5 of 6 pairs (mean difference = 0.065; Wilcoxon signed-rank statistic = 4.0, p = 0.219, reported descriptively — n = 6 is the minimum sample size at which this test was attempted at all). We do not interpret this as evidence of a genuine cell-type-specific dissociation: six paired observations cannot distinguish a real effect from sampling noise.

![Figure 17](../figures/16_pyr_vs_interneuron_paired.png)
**Figure 17. Paired pyramidal-versus-interneuron PV correlation, Analysis D.** PV correlation computed on the identical 6 session-pairs using pyramidal-only versus interneuron-only ensembles, connected by subject/pair. Pyramidal PV correlation exceeds interneuron PV correlation in 5 of 6 pairs (Wilcoxon p = 0.219, descriptive only at this n).

![Figure 18](../figures/17_interneuron_scatter.png)
**Figure 18. NREM duration versus interneuron-ensemble PV correlation, Analysis D.** As Figure 10, using the interneuron-only ensemble tuning curve for the same 6 pairs. Spearman r = −0.429, p = 0.397.

### Summary table

| Analysis | Dataset(s) | Species | n pairs | n subjects | Spearman r | p | Survives RV-coefficient cross-check? |
|---|---|---|---|---|---|---|---|
| A | DANDI:000213 | Rat | 20 | 3 | −0.168 | 0.478 | No |
| B | DANDI:001695 | Mouse | 6 | 3 | −0.829 | 0.042 | — (not independently cross-checked at this n) |
| C | Combined (A+B) | Rat + Mouse | 26 | 6 | −0.246 | 0.225 | No |
| D | DANDI:001695, interneuron ensemble | Mouse | 6 | 3 | −0.429 | 0.397 | — |

All 18 figures (Figures 1–18) are reproduced in full above, inline with the analysis they illustrate; source PNG files are also provided individually in `figures/` for reuse.

---

## 4. Discussion

We set out to test, on chronic multi-day open archival data not collected for this purpose, whether the sleep-dependent hippocampal decorrelation reported by Valdivia et al. (2024) and Bollmann et al. (2025) is also detectable across repeated, many-session exposures to the same environment. Across two independent datasets, two labs, two rodent species, and four converging analyses, the point estimate was consistently in the same direction as the published literature (more NREM associated with lower PV correlation), but no single analysis reached a level of statistical or methodological robustness that we consider adequate to claim a real effect. The combined analysis (Section 3.3) is, by construction, the best-powered estimate available here, and it remains non-significant and fails an independent, cell-identity-free cross-check.

We think three factors most plausibly explain the gap between this result and the published literature, rather than a straightforward failure to replicate: **(1) Sample size.** Even combined, 26 pairs from 6 subjects is small relative to what would be needed to reliably detect an effect of the magnitude reported by Valdivia et al. (2024; r = −0.55) — our own leave-one-subject/dataset-out analyses show the point estimate moving by more than that magnitude depending on which single animal is included. **(2) A genuinely different design.** Both prior studies isolate one learning episode and one, relatively short, intervening sleep opportunity; the present analysis instead uses many repeated exposures to an already-familiar environment, closer to the classic representational-drift literature (Ziv et al., 2013) than to a single-trial memory-consolidation design. Valdivia et al. (2024) explicitly report that their decorrelation effect required the environment to change between encoding and retrieval and was absent in a stable, repeated open-field control — which, if it generalizes, would predict a weaker or absent effect in precisely the repeated-familiar-environment design used here, and may be the single best explanation for our null result. **(3) Methodological compromises forced by archival, not purpose-collected, data.** Most importantly, neither dataset supports cell-identity tracking across days, so our PV correlation is necessarily an ensemble-average proxy rather than the literal, cell-matched population vector correlation used in the primary literature; DANDI:001695's behavioral epochs had to be inferred from position-tracking coverage rather than read from an experimenter-labeled table; and the NREM-between-sessions measure is a lower bound that misses any unmonitored home-cage sleep.

We do not believe these results either confirm or meaningfully disconfirm the sleep-dependent decorrelation hypothesis. We report them because we believe the more valuable contribution of this work is not any single correlation coefficient, but (a) a demonstration that this class of question can be asked of open, chronic, multi-day archival hippocampal data at all, given careful and fully disclosed handling of the gaps such data invariably has relative to a purpose-built experiment, and (b) a fully reusable, extensible pipeline — including the quality-control logic that determines which sessions and subjects are even usable, which we found to be a substantial and non-obvious part of the work (Section 2.2) — that can be pointed at additional sessions as DANDI:001695 and similar dandisets grow, or extended to other open datasets meeting the same criteria (e.g., DANDI:000061, pending pose-extraction from its raw video; DANDI:000629, pending independent spike-sorting of its raw tetrode data, which would offer 27–30 recording days per animal with an explicit sleep–run–sleep–run epoch design).

### 4.1 Relationship to prior work

Both comparator studies used carefully controlled, single-episode designs well suited to isolating a sleep-dependent consolidation effect; the present study trades that precision for scale (multiple animals, multiple labs, many repeated sessions) and for reliance on data collected for other scientific purposes. We see this as complementary rather than competing evidence: a properly powered version of the present design — ideally purpose-collected, with explicit epoch labeling, cross-day cell tracking (e.g., via chronic, depth-stable silicon probes or two-photon imaging), and EMG-supported sleep scoring — would be well positioned to test whether Valdivia et al.'s (2024) finding that decorrelation requires environmental change also holds, or fails to hold, across many repeated sessions in an already-familiar environment.

---

## 5. Limitations

1. **Untracked sleep.** Between-session NREM duration is a lower bound built only from recorded epochs; any home-cage sleep outside the NWB recording window is invisible to this analysis and could not be recovered from either dandiset.
2. **No cross-day single-cell identity tracking.** The ensemble-tuning-curve PV correlation is a documented, disclosed alternative to a literal cell-matched population vector correlation, not equivalent to it.
3. **Small sample size throughout.** The largest analysis (combined, n = 26 pairs) remains a pilot-scale sample by conventional standards; all point estimates should be read with wide, unreported (given the sample size, unstable) confidence intervals in mind.
4. **DANDI:001695 epoch boundaries are inferred, not experimenter-labeled**, from position-tracking coverage, unlike DANDI:000213's explicit epoch table.
5. **DANDI:001695's same-environment criterion is a position-footprint heuristic**, adopted because the session-level metadata was uninformative; it is a reasonable, data-driven proxy but not verified ground truth from the original experimenters.
6. **Electrode/probe depth changes across sessions** in DANDI:000213 mean that some of the day-to-day change in recorded ensemble composition is mechanical (different cells sampled) rather than necessarily reflecting biological reorganization — a confound this design cannot fully separate from genuine representational drift.
7. **Species is confounded with dataset** in the combined analysis (every DANDI:000213 subject is rat; every DANDI:001695 subject is mouse); we treat this explicitly as a covariate/interaction rather than assuming it away, but the available data cannot estimate it reliably.
8. **DANDI:000061 could not contribute to the population-vector analysis** because it lacks extracted position tracking (video only); it is used only descriptively.
9. **Multiple comparisons.** Five distinct robustness checks were run per primary analysis without formal correction; we report the full battery rather than a single selected test, and interpret consistency (or its absence) across the battery qualitatively rather than via a corrected p-value.
10. **This is an observational, not causal, test.** Nothing here establishes that NREM sleep causes any observed (or absent) decorrelation, only whether the two are correlated in these data.

---

## 6. Data and Code Availability

Both primary datasets are openly available without restriction from the DANDI Archive: DANDI:000213 (https://dandiarchive.org/dandiset/000213) and DANDI:001695 (https://dandiarchive.org/dandiset/001695); DANDI:000061 (https://dandiarchive.org/dandiset/000061) was used only descriptively (Section 3.0). All analysis code, the full quality-control/session-selection logic, session-pair tables, figures, and a fully executed, reproducible Jupyter notebook reporting every number in this manuscript are openly available at https://github.com/sisyyphuss/nrem-pv-hippocampus (MIT license). Cached intermediate NWB arrays are not stored in the repository — they are regenerated automatically, by streaming directly from DANDI, when the notebook is re-run — but are otherwise fully reproducible from the two dandisets above.

---

## 7. Author Contributions

S. Yusuf conceived the research question, designed the data selection and quality control protocols, developed the computational pipeline, conducted all statistical analyses, and wrote the manuscript.

## 8. Competing Interests

The author declares no competing interests.

## 9. Acknowledgments

We thank David Tingley, György Buzsáki, Gabrielle Girardeau, Ingrid Inema, Joaquin Gonzalez, and Mihaly Voroslakos for making their chronic hippocampal recordings openly available via the DANDI Archive, without which this analysis would not have been possible.
---

## References

Bollmann, L., Baracskay, P., Stella, F., & Csicsvari, J. (2025). Sleep stages antagonistically modulate reactivation drift. *Neuron*. https://doi.org/10.1016/j.neuron.2025.02.025

Leutgeb, S., Leutgeb, J. K., Barnes, C. A., Moser, E. I., McNaughton, B. L., & Moser, M.-B. (2005). Independent codes for spatial and episodic memory in hippocampal neuronal ensembles. *Science*, 309(5734), 619–623. https://doi.org/10.1126/science.1114037

Robert, P., & Escoufier, Y. (1976). A unifying tool for linear multivariate statistical methods: The RV-coefficient. *Journal of the Royal Statistical Society: Series C (Applied Statistics)*, 25(3), 257–265.

Valdivia, G., Espinosa, N., Lara-Vasquez, A., Caneo, M., Inostroza, M., Born, J., & Fuentealba, P. (2024). Sleep-dependent decorrelation of hippocampal spatial representations. *iScience*, 27(6), 110076. https://doi.org/10.1016/j.isci.2024.110076

Ziv, Y., Burns, L. D., Cocker, E. D., Hamel, E. O., Ghosh, K. K., Kitch, L. J., El Gamal, A., & Schnitzer, M. J. (2013). Long-term dynamics of CA1 hippocampal place codes. *Nature Neuroscience*, 16(3), 264–266. https://doi.org/10.1038/nn.3329

Tingley, D., & Buzsáki, G. (2026). Transformation of a Spatial Map across the Hippocampal-Lateral Septal Circuit (Version draft) [Data set]. DANDI Archive. https://dandiarchive.org/dandiset/000213/draft

Gonzalez, J., & Voroslakos, M. (2026). High-density extracellular recordings (Neuropixels/SiNAPS) of Hippocampal-Cortical dynamics during spatial behavior (Version draft) [Data set]. DANDI Archive. https://dandiarchive.org/dandiset/001695/draft

Girardeau, G., Inema, I., & Buzsáki, G. (2026). Reactivations of emotional memory in the hippocampus-amygdala system during sleep (Version draft) [Data set]. DANDI Archive. https://dandiarchive.org/dandiset/000061/draft

---

## Supplementary Material

- `notebooks/NREM_PV_correlation_analysis.ipynb` — full, executed analysis notebook (Analyses A–D), zero execution errors, all figures embedded in-notebook.
- `results/tables/` — final per-pair data tables underlying every result and figure (`session_pairs.csv`, `session_pairs_001695.csv`, `session_pairs_combined.csv`, `session_pairs_001695_interneuron.csv`).
- `results/figures/` — all 18 figures as standalone, publication-ready PNG files.
- `src/` — full computational pipeline source code (`analysis_lib.py`, `stats_lib.py`, `build_pairs.py`, `make_figures.py`, streaming drivers).
