# Inter-Session NREM Sleep and Hippocampal Population-Vector Decorrelation

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Data: DANDI 000213](https://img.shields.io/badge/DANDI-000213-green.svg)](https://dandiarchive.org/dandiset/000213)
[![Data: DANDI 001695](https://img.shields.io/badge/DANDI-001695-green.svg)](https://dandiarchive.org/dandiset/001695)

This repository provides the computational analysis pipeline, executed Jupyter notebook, and curated statistical tables for the study:

> **Testing Whether Inter-Session NREM Sleep Predicts Hippocampal Population-Vector Decorrelation: An Open, Multi-Dataset Replication Attempt**  
> **Saad Yusuf**  
> *Department of Computer Science, Jamia Millia Islamia, New Delhi, India*  
> Correspondence: `syed.saad7527@gmail.com` | ORCID: `0009-0000-2333-0078`  
> Preprint available in [`preprint/`](preprint/nrem_sleep_pv_decorrelation_manuscript.md).

---

## Scientific Overview

Recent studies have reported that non-rapid-eye-movement (NREM) sleep actively destabilizes or decorrelates hippocampal spatial representations across brief post-learning intervals (Valdivia et al., 2024, *iScience*; Bollmann et al., 2025, *Neuron*). Whether this relationship holds across chronic multi-day recording regimes in familiar environments remains untested.

This project reanalyzes chronic dorsal CA1 recordings across two independent open datasets from the DANDI Archive (rat CA1 silicon probes, [DANDI:000213](https://dandiarchive.org/dandiset/000213); mouse CA1 Neuropixels/SiNAPS, [DANDI:001695](https://dandiarchive.org/dandiset/001695)) to evaluate whether intervening NREM sleep duration predicts population-vector (PV) decorrelation across repeated sessions.

### Primary Findings

1. **Consistent Directional Point Estimate:** Across both rat ($r = -0.168$) and mouse ($r = -0.829$) cohorts, intervening NREM duration correlates negatively with population-vector stability (pooled Spearman $r = -0.246$, $n = 26$ session-pairs, 6 subjects).
2. **Robustness and Sensitivity Profiling:** Under conservative within-subject permutation testing (shuffling NREM duration strictly within animal to preserve marginal distributions), the relationship does not reach statistical significance ($p = 0.224$).
3. **Mixed-Effects Covariate Control:** Linear mixed-effects modeling controlling for intervening wake duration, REM duration, and elapsed calendar days shows no independent NREM effect ($\beta = -0.00048, p = 0.667$).
4. **Cell-Type Control:** Pyramidal ensembles show higher baseline day-to-day spatial stability ($r = 0.109$) compared to interneuron ensembles ($r = 0.044$), ruling out generalized tracking drift as an artifact.

---

## Repository Structure

```
nrem-pv-hippocampus/
├── src/                                  # Core computational and analysis modules
│   ├── analysis_lib.py                   # Rate maps, ensemble tuning curves, PV correlation, RV coefficient
│   ├── stats_lib.py                      # Permutation tests, mixed-effects models, leave-one-out metrics
│   ├── build_pairs.py                    # Session-pair aggregation and intervening sleep accounting
│   ├── make_figures.py                   # Publication-ready figure generation routines
│   ├── cache_session.py                  # DANDI Archive streaming and array caching (DANDI:000213)
│   ├── cache_session_001695.py           # DANDI Archive streaming and array caching (DANDI:001695)
│   ├── cache_session_000061_descriptive.py # Descriptive sleep architecture caching (DANDI:000061)
│   ├── session_list.py                   # Curated inclusion metadata and paths (DANDI:000213)
│   ├── session_list_001695.py            # Curated inclusion metadata and paths (DANDI:001695)
│   ├── probe_session.py                  # Metadata quality-control inspection driver (DANDI:000213)
│   ├── probe_session_001695.py           # Metadata quality-control inspection driver (DANDI:001695)
│   └── run_pipeline.py                   # Master analysis pipeline runner and summary reporter
├── notebooks/                            # Executed Jupyter analysis notebooks
│   └── NREM_PV_correlation_analysis.ipynb # End-to-end reproducible analysis notebook (Analyses A–D)
├── results/                              # Figures and precomputed summary tables
│   ├── figures/                          # Publication-ready figures (Figures 01–17)
│   └── tables/                           # CSV session-pair tables and probe inspection logs
├── preprint/                             # Manuscript writeup
│   └── nrem_sleep_pv_decorrelation_manuscript.md
├── requirements.txt                      # Python library dependencies
├── LICENSE                               # MIT License
└── README.md
```

---

## Description of Analysis Modules (`src/`)

- **`analysis_lib.py`:** Computes occupancy-normalized 2D rate maps, extracts ensemble tuning curves across valid spatial bins, and calculates population-vector Pearson correlations and Robert-Escoufier RV coefficients.
- **`stats_lib.py`:** Implements pooled Spearman/Pearson rank tests, linear mixed models (`statsmodels`), leave-one-subject-out (LOSO) cross-validation, and within-subject permutation null models (5,000 permutations).
- **`build_pairs.py`:** Iterates through consecutive same-environment sessions, aligns spatial coordinates, quantifies intervening NREM/REM/Wake durations, and constructs pair-level analysis tables.
- **`run_pipeline.py`:** Command-line driver that executes the complete statistical battery across single-dataset and pooled cohorts and prints formal summary tables.
- **`make_figures.py`:** Standalone plotting utilities using an Okabe-Ito colorblind-safe palette for all manuscript figures.

---

## Data Sources & Inclusion Criteria

- **[DANDI:000213](https://dandiarchive.org/dandiset/000213):** Tingley et al. (Buzsáki Lab). Rat dorsal CA1 and lateral septum silicon-probe recordings. 3 subjects retained after QC (`DT2`, `DT5`, `DT7`), yielding 20 same-task session pairs.
- **[DANDI:001695](https://dandiarchive.org/dandiset/001695):** Gonzalez et al. Mouse CA1 Neuropixels/SiNAPS recordings. 3 subjects retained (`M01`, `M02`, `M03`), yielding 6 same-configuration session pairs verified via spatial tracking footprint clustering.
- **[DANDI:000061](https://dandiarchive.org/dandiset/000061):** Girardeau et al. (Buzsáki Lab). Rat CA1/BLA recordings; evaluated descriptively for sleep-state architecture.

---

## Reproduction & Usage

### 1. Environment Setup

```bash
git clone https://github.com/sisyyphuss/nrem-pv-hippocampus.git
cd nrem-pv-hippocampus
pip install -r requirements.txt
```

### 2. Running the Pipeline

Execute the master analysis driver to reproduce all statistical tests from cached session tables:

```bash
python src/run_pipeline.py
```

### 3. Interactive Notebook

Launch Jupyter to inspect step-by-step analyses, diagnostic checks, and figure plots:

```bash
jupyter lab notebooks/NREM_PV_correlation_analysis.ipynb
```

---

## Citation

If you use this pipeline or data tables in your research, please cite:

```bibtex
@article{yusuf2026nrem,
  title={Testing Whether Inter-Session NREM Sleep Predicts Hippocampal Population-Vector Decorrelation: An Open, Multi-Dataset Replication Attempt},
  author={Yusuf, Saad},
  journal={bioRxiv},
  year={2026},
  note={Preprint}
}
```

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
