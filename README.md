# MPNc model: Application reproduction and soft-closure simulation

This repository reproduces the **Application** (Section 7) of the paper
*Identifiability and Observational Equivalence in Soft-Closure Multinomial
N-Mixture Models* on the Chestnut-sided Warbler data (`cswa`), and provides a
simulation CLI (Section 6) for the four models MPN0, MPN1, MPNc and CRTE.

The whole analysis is written in Python; R is not required. The source `.rda`
file is read with the `rdata` package when it is installed, and otherwise with
the dependency-free decoder `rda_reader.py`; the tests confirm that both give
the same data. The detailed report is in
[`APPLICATION_REPORT.md`](APPLICATION_REPORT.md).

## Status of the included results

The original results were produced with Python 3.12.14 (NumPy 2.3.5,
SciPy 1.17.0) and rerun with Python 3.13.16 (NumPy 2.5.3, SciPy 1.18.1). The
results in `results/` come from the most recent rerun, with Python 3.12.3
(NumPy 2.2.6, SciPy 1.16.3) and the `rda_reader` decoder
(`results/run_manifest.json`), and agree with the earlier runs to at least 6
decimal places. All 8 fits succeed; **77/77** numbers from Section 7 of the
manuscript (estimates, logLik, AIC, weights, confidence intervals, third-order
diagnostic, composition test, total area) match after rounding
(`results/paper_comparison.csv`).

| Analysis | Best model by AIC | AIC | Abundance estimate | `r` per minute | Dependence parameter |
| --- | --- | ---: | ---: | ---: | ---: |
| Area offset | MPNc | 353.666008 | 1.083596 individuals/ha | 0.373167 | `c = 0.456632` |
| Intercept only | MPNc | 356.120210 | 0.589969 individuals/point/visit | 0.371773 | `c = 0.456554` |

These are results on the real data, not a simulation smoke test. See
`results/verification_summary.json` for convergence checks and
`results/paper_comparison.csv` for each individual comparison.

## Installation

The library versions used for the included run are pinned in
`requirements.txt`. `rdata` and `xarray` are optional: without them the code
falls back to `rda_reader.py`.

```bash
git clone https://github.com/huynhhuudinh/MPNc.git
cd MPNc
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

On Windows, activate the environment with `.venv\Scripts\activate`. The scripts
resolve their default paths from their own location, so they can also be called
from another directory using an absolute path.

## Rerunning the Application

The full workflow runs the tests, prepares the data, fits the models and builds
the report:

```bash
python -m unittest -v test_application.py
python run_application.py --starts 8
python build_report.py
```

`run_application.py` runs both specifications (`area` and `intercept`) for all
four models. Available options:

```text
--data-dir PATH   data directory; default data/
--out PATH        results directory; default results/
--download        re-download the three pinned source files and check SHA-256
--prepare-only    only check and export the standardized data; do not fit models
--starts N        number of starting points; default 8, must be positive
```

To keep the included results untouched:

```bash
python run_application.py --out results_new --starts 8
python build_report.py --results results_new
```

A new run only overwrites output files with the same names in the chosen
directory; it does not delete other files there. `build_report.py --results ...`
reads results from the given directory, writes the check tables into that same
directory, and always updates `APPLICATION_REPORT.md` in the project directory.

### Data

`data/raw/` contains `cswa.rda`, its documentation `data_cswa.Rd`, the package
metadata `DESCRIPTION` and `provenance.json`. The source is pinned to a commit
of the CRAN mirror of `AHMbook`; SHA-256 hashes are checked before every read.
No network access is needed when using the included files.

```bash
# Re-download the pinned source, check the hashes, then run the analysis
python run_application.py --download --starts 8

# Only verify the source and regenerate data/*.csv, data/*.npz, data_audit.json
python run_application.py --prepare-only
```

The data have 43 points, 3 visits and removal intervals of 2, 3 and 5 minutes.
There are 128 observed visits; visit 1 of point `darty` is missing entirely and
is excluded from the likelihood. In `counts_long.csv` the row for the missing
visit is kept with an empty `count` and `observed=False`; in `analysis_data.npz`
the value stays `NaN` and is not replaced by 0. The area covariate is joined by
`point_id`, not by row position.

Standardized data files:

| File | Contents |
| --- | --- |
| `data/counts_long.csv` | Counts in long format, with observation flag and plot area |
| `data/site_covariates.csv` | Per-point covariates, aligned by identifier |
| `data/spot_mapping.csv` | Spot-mapping data by patch |
| `data/analysis_data.npz` | Arrays `y`, `area_ha` and `point_id` used in the analysis |
| `data/data_audit.json` | Data axes, missingness, joins, units and check totals |

### Models and units

- **MPN0**: abundance independent across visits; this is a working likelihood
  when the real data are dependent. Its primary standard errors use the
  point-clustered sandwich.
- **MPN1**: abundance fully closed across visits.
- **MPNc**: a shared Poisson component with mean `c * mu_i` and a visit-specific
  component with mean `(1 - c) * mu_i`.
- **CRTE**: a Poisson superpopulation with mean `mu_i / phi`; availability is
  independent across visits with probability `phi`. `mu_i` is the mean number
  present at a single visit.

In the `area` analysis, `mu_i = rho * plotArea_i` and `rho` is in
individuals/ha. In the `intercept` analysis, `mu_i = mu` and the unit is
individuals/point/visit. `r` is always the per-minute detection probability.
The `dep` column is `c` for MPNc and `phi` for CRTE, and is fixed at 0 and 1 for
MPN0 and MPN1 respectively.

Likelihoods are computed on the log scale and keep the factorial/multinomial
constants so that absolute logLik and AIC values are reproduced. Sums over latent
abundance use a tail threshold plus a stricter check; the numerical limits are
set in `FitConfig` in `softclosure.py`.

### Application outputs

| File | Contents |
| --- | --- |
| `results/model_results.csv` | Estimates, logLik, AIC, weights, SEs, CIs and diagnostics for the 8 fits |
| `results/area_*.json` | Per-model details for the area-offset specification |
| `results/intercept_*.json` | Per-model details for the intercept-only specification |
| `results/site_loglik.csv` | Log-likelihood contribution by point and model |
| `results/third_order.json` | Area-based third-order residual diagnostic and pooled version for audit |
| `results/composition_gof.json` | Goodness-of-fit test of the removal composition (X², G², bootstrap p) |
| `results/run_manifest.json` | Software versions, configuration, data source and hashes of the code that ran |
| `results/paper_comparison.csv` | Comparison of 77 Section 7 numbers with rounding tolerances |
| `results/confidence_intervals.csv` | Log/logit CIs and natural-scale Wald intervals, with covariance type |
| `results/verification_summary.json` | Summary of convergence, numerical stability and comparisons |
| `results/tests.log`, `results/run.log` | Logs of the included test run and analysis run |

In the JSON/CSV files, the `mu_` prefix from the shared numerical module refers
to the first abundance parameter: read it as `rho` for the `area` analysis and
`mu` for `intercept`. `confidence_intervals.csv` already renames parameters by
estimand. The `primary_se` field indicates which covariance is reported as the
primary result.

## Running the simulation

`softclosure.py` both provides the numerical functions reused by
`run_application.py` and is a standalone simulation CLI for designs with equal
removal intervals. Available suites are `smoke` (9 scenarios), `core` (24),
`reverse` (24), `small` (6), `stress` (6) and `full` (60).

The simulation summary produces PNG/PDF figures and additionally needs
`matplotlib` (not in `requirements.txt`, which covers the Application workflow):

```bash
python -m pip install matplotlib
python softclosure.py run --suite smoke --reps 10 --workers 2 --out simulation_smoke
```

Without `--reps`, the default is 10 for `smoke` and 500 for the other suites.
`--seed`, `--starts`, `--maxiter` and `--no-intervals` are also available;
`--workers` must be between 1 and 8. Each replicate is written separately under
`OUT/replicates/`, so an interrupted run can be resumed as long as the
fingerprint of the code, seed, scenario and configuration still matches. Do not
reuse the same output directory after changing any of these.

Rebuild the summary tables and figures from the completed replicates, then build
the tables and figures used in the paper (fit-rate table, model-selection
tables, misspecification cell, and the two figures `full_misspecification` and
`full_model_selection`) into `OUT/paper/`:

```bash
python softclosure.py report --out simulation_smoke
python paper_simulation_outputs.py --out simulation_smoke
```

The paper's full run: `python softclosure.py run --suite full --out simulation_full`
(default seed 20260910, 500 replicates, 3 starting points). It takes about
1.5–4 seconds per data set on one core, i.e. a few hours with 8 workers.

Simulation outputs are `manifest.json`, `fits.csv`, `parameter_summary.csv`,
`model_selection.csv`, `selection_summary.csv`, `third_order.csv`,
`third_order_summary.csv`, `run_summary.json` and the `rmse_*` figures in PNG/PDF
format. The smoke suite only checks that the pipeline runs; it is not full
simulation evidence.

The included `simulation_full/` holds the aggregated outputs of the paper's full
run. The 30,000 per-replicate files in `simulation_full/replicates/` (about 186 MB) are
not included in this repository; their contents are aggregated in `fits.csv`,
`model_selection.csv` and `third_order.csv`, and they can be regenerated with the
command above.

## Code layout

| File | Role |
| --- | --- |
| `run_application.py` | Source verification, data standardization, the 8 model fits, third-order diagnostic, composition test |
| `rda_reader.py` | Dependency-free `.rda` decoder (RDX2/RDX3, XDR) |
| `softclosure.py` | Shared likelihood/numerical core and resumable simulation CLI |
| `paper_simulation_outputs.py` | Builds the paper's simulation tables and figures from a results directory |
| `test_application.py` | 10 tests: removal probabilities, data axes, missingness, the four likelihoods against the brute-force multinomial sum, the independent decoder, the composition test, the working profile |
| `build_report.py` | Checks 77 numbers against the paper, builds the CI table, verification summary and Markdown report |
| `PACKAGE_MANIFEST.json` | Snapshot of sizes and SHA-256 hashes of the delivered files, including `simulation_full/` (except the 30,000 files in `simulation_full/replicates/`, for which only the count and total size are recorded) |
| `LICENSE` | GPL-3.0-or-later license |

## Scope and interpretation

The Application workflow only reproduces the intercept and area specifications;
it does not add covariates for time, date, observer, habitat, patch effects or
overdispersion. The data have no individual identities, so `c` does not directly
show that the same individuals return. As in the paper, the analysis treats the
43 points as independent even though they belong to 38 patches; this assumption
is not tested here. AIC weights are not probabilities that a biological mechanism
is true.

The data source credits David King and colleagues; the source documentation
cites Chandler, King & Chandler (2009), Chandler, Royle & King (2011), and Kéry &
Royle (2021). The AHMbook package is licensed GPL (>= 3); redistribution must
keep the files and provenance in `data/raw/`. The project's Python code is
provided under GPL-3.0-or-later; see `LICENSE`.
