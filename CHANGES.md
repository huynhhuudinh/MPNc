# Changes relative to the original package (2026-10-12)

## Review findings

No mathematical errors were found in the likelihoods, fitting algorithm,
confidence intervals, third-order diagnostic or simulation summaries. Checks
performed:

- The four likelihoods were compared with direct summation over every original
  multinomial allocation (existing test), and with the bivariate Poisson
  distribution in the two-visit case.
- The `working_fit` profile solution coincides with the numerical maximum of the
  independence likelihood (log-likelihood difference < 1e-9).
- The data generator reproduces the theoretical means, covariances and
  third-order cross-cumulants of both MPNc and CRTE (400,000 sites).
- The standardized data were decoded again independently from `cswa.rda`
  (without `rdata`) and match exactly, including the name-based join of the last
  three points.
- Application rerun: agrees with the original run to ≥ 6 digits; 77/77 numbers
  from Section 7 of the manuscript match after rounding.
- The simulation was rerun for exactly the three published cells with the same
  seed 20260910 (500 replicates per cell): Table 3, the correct-selection rate of
  58.3% (MCSE 0.0220), the median ΔAIC of 0.172 and R̂₃ (0.183; SD 0.300) match
  exactly; the results for the two cells of Table 5 are in `verification/`.

## Changes

1. **`rda_reader.py` (new).** A dependency-free RDX2/RDX3 decoder.
   `run_application.read_source` uses `rdata` when available and otherwise this
   decoder; `run_manifest.json` records which decoder was used. Reason:
   `rdata`/`xarray` cannot always be installed, and the data-reading step needs
   an independent check.
2. **Composition test (`composition_gof`, new).** The revised manuscript reports
   X², G² and p for the pooled composition (Section 7.3), but the original code
   did not compute them. The new function writes `results/composition_gof.json`,
   including a parametric bootstrap p (4,999 draws, seed 20261012) because the
   smallest expected cell is only 6.6.
3. **`build_report.py` rewritten.**
   - Removed the outdated "Points to fix: CI" section (the manuscript now uses
     log/logit CIs).
   - Every statement in the report is computed from the results instead of being
     hard-coded (areas, re-joined points, direction of the sandwich/information
     SE comparison, model with the smallest AIC).
   - The comparison was extended from 22 to 77 numbers, each with a tolerance
     based on the number of printed digits; mismatches are reported explicitly.
4. **`paper_simulation_outputs.py` (new).** The Data Availability section states
   that the simulation tables and figures are regenerated from the summaries, but
   the original package had no script for the two figures
   `full_misspecification`, `full_model_selection` and the pooled tables. The new
   script builds them from a `softclosure.py` results directory.
5. **Tests: 7 → 10.** Added: the independent decoder against the standardized
   data; the composition test (value on the real data and a near-perfect-fit
   case); the working profile coinciding with the numerical MLE.
6. **README.** Fixed the title to match the manuscript, removed local machine
   paths, described the new files, how to build the paper's tables/figures and
   the estimated time for the full run.

`softclosure.py`, `rda_reader.py` and `run_application.py` do not change any
likelihood formula or optimizer setting; the numerical results are unchanged.

## Translation to English (2026-11-26)

`README.md`, `CHANGES.md`, `APPLICATION_REPORT.md` and the report template in
`build_report.py` were translated from Vietnamese to English. Only the report
text changed; the computations in `build_report.py` are unchanged, and
regenerating the report leaves every file in `results/` byte-identical.
`results/run_manifest.json` still records the SHA-256 of `build_report.py` as it
was at the time of the run, so that hash no longer matches the current file.
