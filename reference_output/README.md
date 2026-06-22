# Reference output

These are the exhibits as they appear in the paper, provided so that a replicator
can compare their own run against the published results:

| File | Paper exhibit | Produced by |
| --- | --- | --- |
| `Table 1.csv` | Table 1 (main simulation, non-stationary DGP) | `simulation/Main Simulation Results (Table1&2)/funs/simulation_summary_table1&2.R` |
| `Table 2.csv` | Table 2 (main simulation, stationary DGP) | `simulation/Main Simulation Results (Table1&2)/funs/simulation_summary_table1&2.R` |
| `Table 3.csv` | Table 3 (bandwidth choice) | `simulation/Bandwidth Choice (Table 3)/funs/simulation_summary_table3.R` |
| `Table 4.xlsx` | Table 4 (Sequeira 2016 application) | `application/funs/summary_table.R` |
| `Figure 1.pdf` | Figure 1 (power) | `simulation/Power Analysis (Figure 1)/funs/simulation_summary_power.R` |

How to use: after running a scenario's driver(s) and its `simulation_summary_*`
script (or the application scripts), compare the file written to that scenario's
`tables/` or `figures/` directory against the corresponding file here.

Note on numerical agreement: the local-polynomial routines (pseudo)invert local
moment matrices that can be near-singular, so the simulation results may differ
slightly across operating systems, compilers, and BLAS/LAPACK builds — see the
"Reproducibility" section of the top-level `README.md`. The CSV tables here are
the exact reference; small last-digit differences in a fresh run are expected and
do not change the qualitative conclusions.
