# Difference-in-Differences with Compositional Changes — Replication Package

Replication code for *Difference-in-Differences with Compositional Changes*.

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
  bribe payments at the port of Maputo (Sequeira, 2016, *AER*), Table 4.

---

## Repository structure

```
.
├── README.md
├── LICENSE
├── install_packages.R
├── .gitignore
│
├── simulation/                              # Section 5 — Monte Carlo study
│   ├── core/                                # shared backend, sourced by every driver
│   │   ├── np_lp.cpp, np_lp_rcv.cpp         #   local-polynomial routines (C++ / Rcpp)
│   │   ├── dgps_drdid.R                     #   data-generating process; se.bound()
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
    ├── funs/  main_sequeira.R, main_sequeira_twfe.R, sequeira_table.R,
    │          pre_process.R, locpol_funs.R, dr_did.R, brackets.R, np_lp.cpp
    ├── data/   Bribes_Regression.dta        # obtain separately (git-ignored)
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

## Software requirements

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

---

## Setup — read before running

1. **Set the local path.** Each driver and summary script sets an absolute path
   near the top (`address`, `address0`, or `base_dir`) flagged with
   `# EDIT to your local path`. Change it to point to your copy of the relevant
   folder before running.
2. **Run the application from `application/`.** The application backend is
   sourced with a relative path (`funs/np_lp.cpp`), so the working directory
   must be the `application/` folder when you run the Sequeira scripts.
3. **Provide the data.** The application data file is not included — see
   *Data availability* below — place `Bribes_Regression.dta` in
   `application/data/`.

---

## Reproducing the simulations (Section 5)

For every scenario the workflow is the same: run the **`main_*`** driver(s) to
generate per-job results into `results/`, then run the matching
**`simulation_summary_*`** script to aggregate them into the table/figure.

> **Runtime.** A full run is very long — on the order of three weeks of CPU time.
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

With the working directory set to `application/`:

1. `funs/main_sequeira.R` — doubly robust DiD estimates (Panel B) and the
   Hausman-type stationarity test (Panel C).
2. `funs/main_sequeira_twfe.R` — two-way fixed-effects estimates (Panel A).
3. `funs/sequeira_table.R` — assembles Table 4 from the saved results
   (`results/table_ml_cluster.xlsx`).

Inference uses a cluster bootstrap (clustering on HS 4-digit codes); the seed is
fixed in each script.

---

## Outputs and version control

`results/`, `tables/`, and `figures/` hold generated output and are **not**
tracked (see `.gitignore`); only the empty directory skeleton is kept via
`.gitkeep` so the scripts have somewhere to write after a fresh clone.
Re-running the scripts regenerates everything.

## Reproducibility

Each script fixes its random seed (`seed1`), and the reported results were
produced with R 4.3.3 and the package versions in `install_packages.R` (notably
`np` 0.60.17). Several components are estimated nonparametrically — local-
polynomial generalized propensity scores and outcome regressions, with
bandwidths chosen by cross-validation or plug-in rules — through routines
compiled from C++ via `RcppArmadillo`.

These differences arise in the local-polynomial step: at each evaluation point
it (pseudo-)inverts a local weighted moment matrix (`X' W X`), which can be
near-singular when few observations carry appreciable weight there. Inverting an
ill-conditioned matrix amplifies small floating-point differences across
operating systems, compilers, and BLAS/LAPACK builds — a numerical feature of
local estimation, that leaves the qualitative conclusions unchanged.

This is most visible in the plug-in column of Table 3: it uses a frequency-based
kernel (no smoothing) for the discrete covariates, so each fit conditions on an
exact discrete cell, giving smaller local samples and more often near-singular
matrices than the cross-validated (LOOCV/RCV) variants, which smooth across
categories. Matching the reported R and package versions (notably `np`) and
toolchain gives the closest agreement; the shipped result files are the exact
reference.

### Troubleshooting

* **`pinv(): svd failed`.** On some platforms / linear-algebra backends the
  pseudo-inverse in the local-polynomial routines can fail to converge with the
  standard SVD method. If you hit this error, open `np_lp.cpp` (in
  `simulation/core/` and `application/funs/`) and change the affected `pinv()`
  calls from `pinv(A, DBL_SMALL, "std")` to `pinv(A, DBL_SMALL, "dc")`, then recompile. This
  falls back to Armadillo's default pseudo-inverse and lets the routine run on
  the affected machines. (On machines where the original form already runs, no
  change is needed.)

## Data availability

The application uses the replication data from Sequeira (2016), *Corruption,
Trade Costs, and Gains from Tariff Liberalization*, American Economic Review.
Obtain `Bribes_Regression.dta` from the journal's replication archive and place
it in `application/data/`. The data are not redistributed in this repository.

## License

Code is released under the MIT License (see `LICENSE`). The Sequeira (2016) data
remain subject to their original license/terms.

## Citation

If you use this code, please cite the paper:

> [Pedro H.C. Sant’Anna, Qi Xu] (2026). *Difference-in-Differences with Compositional Changes.*
> Journal of Econometrics, Volume 253, 2026, 106147, ISSN 0304-4076, https://doi.org/10.1016/j.jeconom.2025.106147.

