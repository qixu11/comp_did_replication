# Difference-in-Differences with Compositional Changes — Replication Package

This repository contains the replication code for [Sant'Anna and Xu (2026)](https://doi.org/10.1016/j.jeconom.2025.106147), *"Difference-in-Differences with Compositional Changes"*.

The package implements doubly robust (DR) DiD estimators for repeated cross
sections that remain valid when the covariate distribution shifts across
periods/groups ("compositional changes" / non-stationarity). Nuisance
functions — the generalized propensity score (GPS) and the outcome regression
(OR) — are estimated by local-polynomial regression, and a Hausman-type test
is used to test for covariate stationarity.

It reproduces:

* **Section 5 — Monte Carlo study:** Tables 1 & 2 (main results), Table 3
  (bandwidth choice), and Figure 1 (power).
* **Section 6 — Empirical application:** the effect of tariff reductions on
  bribe payments at the port of Maputo in Mozambique ([Sequeira, 2016, *AER*](http://dx.doi.org/10.1257/aer.20150313)), Table 4.

---

## Repository structure

```text
.
├── README.md
├── LICENSE                                  # MIT — covers the code in this repository
├── install_packages.R
├── .gitignore
│
├── reference_output/                        # the paper's exhibits, for comparison
│   ├── Table 1.csv, Table 2.csv, Table 3.csv, Table 4.xlsx, Figure 1.pdf
│   └── README.md
│
├── simulation/                              # Section 5 — Monte Carlo study
│   ├── core/                                # shared backend, sourced by every driver
│   │   ├── np_lp.cpp, np_lp_rcv.cpp         #   local-polynomial routines (C++ / Rcpp)
│   │   ├── dgps_drdid.R                     #   data-generating process
│   │   ├── locpol_funs.R                    #   GPS/OR local-poly fitting + CV bandwidths
│   │   ├── dr_did.R                         #   DR DiD + TWFE estimators + bootstrap
│   │   └── bw_plugin.R                      #   plug-in bandwidth selection
│   │
│   ├── Main Simulation Results (Table1&2)/
│   │   ├── funs/  main_drdid.R, sim_drdid.R, simulation_summary_table1&2.R
│   │   ├── results/                         # generated output (git-ignored)
│   │   └── tables/                          # generated output (git-ignored)
│   │
│   ├── Bandwidth Choice (Table 3)/
│   │   ├── funs/  main_bdw_{loocv,rcv,plugin}.R, sim_bdw_{loocv,rcv,plugin}.R,
│   │   │          simulation_summary_table3.R
│   │   ├── results/{loocv,rcv,plugin}/      # generated output (git-ignored)
│   │   └── tables/                          # generated output (git-ignored)
│   │
│   └── Power Analysis (Figure 1)/
│       ├── funs/  main_drdid_power.R, sim_power.R, simulation_summary_power.R
│       ├── results/                         # generated output (git-ignored)
│       └── figures/                         # generated output (git-ignored)
│
└── application/                             # Section 6 — Sequeira (2016) application
    ├── funs/  main_drdid.R, main_twfe.R, summary_table.R,
    │          pre_process.R, locpol_funs.R, dr_did.R, brackets.R, np_lp.cpp
    ├── data/  Bribes_Regression.dta         # Sequeira (2016); bundled under CC BY 4.0
    │          LICENSE-Sequeira-2016-AEA.txt #   data license + attribution
    ├── results/                             # generated output (git-ignored)
    └── tables/                              # generated output (git-ignored)
```

Within each simulation scenario the `funs/` folder follows the same pattern:
a **`main_*`** driver (sets parameters and loops the Monte Carlo), a **`sim_*`**
inner loop (one replication's estimation, sourced by the driver), and a
**`simulation_summary_*`** script that aggregates the saved results into the
paper's table/figure. The simulation drivers all source the shared backend in
`simulation/core/`; the application carries its own copy of the backend in
`application/funs/`.

---

## List of exhibits

Each paper exhibit is produced by a driver (or set of drivers) followed by a
summary script. Reference copies of every output are in `reference_output/`.

| Exhibit | Section | Driver(s) | Summary script | Output file |
| --- | --- | --- | --- | --- |
| Tables 1 & 2 | 5 | `main_drdid.R` | `simulation_summary_table1&2.R` | `tables/Table 1.csv`, `Table 2.csv` |
| Table 3 | 5 | `main_bdw_loocv.R`, `main_bdw_rcv.R`, `main_bdw_plugin.R` | `simulation_summary_table3.R` | `tables/Table 3.csv` |
| Figure 1 | 5 | `main_drdid_power.R` | `simulation_summary_power.R` | `figures/Figure 1.pdf` |
| Table 4 | 6 | `main_drdid.R`, `main_twfe.R` (in `application/`) | `summary_table.R` | `tables/Table 4.xlsx` |

---

## Software and computational requirements

* **R** (tested on 4.3.3) with a working C++ toolchain, because the
  local-polynomial backend is compiled at run time via `Rcpp::sourceCpp()`:
  * macOS — Xcode Command Line Tools (`xcode-select --install`)
  * Windows — Rtools
  * Linux — a system C++ compiler (e.g. `g++`)
* **R packages:** `Rcpp`, `RcppArmadillo`, `roptim`, `foreach`, `doParallel`,
  `doRNG`, `np`, `nnet`, `MASS`, `glmnet`, `haven`, `dplyr`, `ggplot2`,
  `reshape2`, `patchwork`, `scales`, `openxlsx`.
  (`MASS` and `nnet` ship with R.)

Install everything with:

```r
source("install_packages.R")
```

* **Operating system:** platform-independent. The code has no OS-specific
  dependencies and runs on macOS, Windows, and Linux given the toolchain above.
* **Hardware and runtime:** the empirical application (Table 4) and every
  `simulation_summary_*` aggregation script run in minutes on a standard laptop
  (≈ 8 GB RAM). The full Monte Carlo study is heavy — on the order of several weeks
  of single-core CPU time in total — and was originally run as parallel cluster
  jobs; see *Reproducing the simulations* for how to scale it down for a quick
  check.
* **Disk:** the bundled data and the reference outputs are < 1 MB combined; a
  full simulation run writes a few hundred MB of intermediate files to `results/`.

---

## Setup — read before running

1. **Run everything from the repository root.** Every driver and summary script
   uses repository-root-relative default paths (`address`, `address0`, `base_dir` /
   `app_dir`), so with the working directory set to the repo root they resolve
   automatically — no editing required. (Advanced: you may instead set those
   variables near the top of each script to absolute paths, in which case the
   working directory no longer matters.)
2. **Data is included.** The application dataset (`Bribes_Regression.dta`) is
   bundled in `application/data/` under its original license — no download is
   required. See *Data availability and provenance* below for source and terms.

---

## Reproducing the simulations (Section 5)

For every scenario the workflow is the same: run the **`main_*`** driver(s) to
generate per-job results into `results/`, then run the matching
**`simulation_summary_*`** script to aggregate them into the table/figure.

> **Runtime.** A full run is very long — on the order of several weeks of CPU time.
> Each driver loops over jobs `0`–`9`; in practice each job is run as a separate
> (parallel) cluster task. For a quick check, reduce `nrep` and/or the job range
> at the top of the driver.

### Tables 1 & 2 — main simulation

1. Run `simulation/Main Simulation Results (Table1&2)/funs/main_drdid.R`
2. Run `simulation_summary_table1&2.R`
   → `tables/Table 1.csv`, `Table 2.csv`

### Table 3 — bandwidth choice

1. Run all three drivers: `main_bdw_loocv.R`, `main_bdw_rcv.R`,
   `main_bdw_plugin.R` (LOOCV, RCV, and plug-in bandwidth selectors)
2. Run `simulation_summary_table3.R` → `tables/Table 3.csv`

### Figure 1 — power

1. Run `simulation/Power Analysis (Figure 1)/funs/main_drdid_power.R`
2. Run `simulation_summary_power.R` → `figures/Figure 1.pdf`

---

## Reproducing the application (Section 6, Table 4)

With the working directory set to the repository root, and
`application/data/Bribes_Regression.dta` in place:

1. `application/funs/main_drdid.R` — doubly robust DiD estimates (Panel B) and the
   Hausman-type stationarity test (Panel C). Writes
   `application/results/result_sequeira_drdid.RData`.
2. `application/funs/main_twfe.R` — two-way fixed-effects estimates (Panel A).
   Writes `application/results/result_sequeira_twfe.RData`.
3. `application/funs/summary_table.R` — assembles Table 4 from those two saved
   `.RData` files into `application/tables/Table 4.xlsx`.

Inference uses a cluster bootstrap (clustering on HS 4-digit codes); the seed is
fixed in each script.

---

## Outputs and version control

`results/`, `tables/`, and `figures/` hold generated output and are **not**
tracked (see `.gitignore`); only the empty directory skeleton is kept via
`.gitkeep` so the scripts have somewhere to write after a fresh clone.
Re-running the scripts regenerates everything. The paper's exhibits are kept
(tracked) under `reference_output/` so a replicator can compare a fresh run
against the published numbers.

## Reproducibility

Each script fixes its random seed (`seed1`), and the reported results were
produced with R 4.3.3 and, in particular, `np` 0.60.17. `install_packages.R`
installs the current CRAN build of each package rather than pinning exact
versions, so for the closest agreement install the versions noted here. Several components are estimated nonparametrically — local-
polynomial generalized propensity scores and outcome regressions, with
bandwidths chosen by cross-validation or plug-in rules — through routines
compiled from C++ via `RcppArmadillo`.

Small numerical differences arise in the local-polynomial step. At each evaluation point, the procedure computes a (pseudo-)inverse of a local weighted moment matrix, (`X' W X`), which can become nearly singular when few observations receive appreciable weight. Inverting an ill-conditioned matrix amplifies small floating-point differences across operating systems, compilers, and BLAS/LAPACK implementations. This is a numerical feature of local estimation and does not affect the paper's qualitative conclusions.

These differences are most apparent in the plug-in rows of Table 3. The plug-in estimator uses a frequency-based kernel (that is, no smoothing) for the discrete covariates, so each fit conditions on an exact discrete cell. This results in smaller effective local samples and more frequent near-singular matrices than the cross-validated (LOOCV/RCV) variants, which smooth across categories. Matching the reported R and package versions (especially `np`) and using a similar toolchain yields the closest agreement.

### Troubleshooting

* **`pinv(): svd failed`.** On some platforms / linear-algebra backends the
  pseudo-inverse in the local-polynomial routines can fail to converge with the
  standard SVD method. If you hit this error, open `np_lp.cpp` (in
  `simulation/core/` and `application/funs/`) and change the affected `pinv()`
  calls from `pinv(A, DBL_SMALL, "std")` to `pinv(A, DBL_SMALL, "dc")`, then recompile. This
  falls back to Armadillo's default pseudo-inverse and lets the routine run on
  the affected machines. (On machines where the original form already runs, no
  change is needed.)

## Data availability and provenance

The empirical application uses a single dataset, which is **included** in this
repository.

**Rights and license.** The data are © 2016 American Economic Association and are
made available under the Creative Commons Attribution 4.0 International license
(**CC BY 4.0**). They are redistributed here **unmodified** under that license;
the full license text and attribution are in
[`application/data/LICENSE-Sequeira-2016-AEA.txt`](application/data/LICENSE-Sequeira-2016-AEA.txt).
This package and its authors are not affiliated with, nor endorsed by, the
original author or the American Economic Association.

**Dataset list.**

| Dataset | File | Source | Provided here? | License |
| --- | --- | --- | --- | --- |
| Bribe-payment regression sample (Southern African ports) | `application/data/Bribes_Regression.dta` (662 KB) | Sequeira (2016) replication package (openICPSR) | Yes — bundled | CC BY 4.0 |

**Data citations.** The dataset is derived from, and should be cited as:

* Sequeira, Sandra. 2016. "Corruption, Trade Costs, and Gains from Tariff
  Liberalization: Evidence from Southern Africa." *American Economic Review*
  106 (10): 3029–63. <https://doi.org/10.1257/aer.20150313>
* Sequeira, Sandra. 2016. "Replication data for: Corruption, Trade Costs, and
  Gains from Tariff Liberalization: Evidence from Southern Africa." American
  Economic Association [publisher], Inter-university Consortium for Political and
  Social Research [distributor]. <https://doi.org/10.3886/E113044V1>

## License

* **Code** in this repository is released under the MIT License (see `LICENSE`).
* **Data** (`application/data/Bribes_Regression.dta`) is © 2016 American Economic
  Association and is redistributed under CC BY 4.0; see
  [`application/data/LICENSE-Sequeira-2016-AEA.txt`](application/data/LICENSE-Sequeira-2016-AEA.txt)
  and *Data availability and provenance* above.

## Citation

If you use this code, please cite the paper:

> Pedro H.C. Sant’Anna and Qi Xu (2026). *Difference-in-Differences with Compositional Changes.*
> Journal of Econometrics, Volume 253, 2026, 106147, ISSN 0304-4076, https://doi.org/10.1016/j.jeconom.2025.106147.

