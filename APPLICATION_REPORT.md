# Application reproduction results

Ran **8 fits: 4 models × 2 abundance specifications** on the real data. 8/8 fits succeeded. Checked 77 numbers printed in Section 7 of the manuscript (estimates, logLik, AIC, weights, confidence intervals, third-order diagnostic, composition test, total area): 77/77 match after rounding. All match. MPNc has the smallest AIC in the primary analysis; MPNc in the intercept-only analysis.

## Data

Source: [AHMbook cswa.rda](https://github.com/cran/AHMbook/blob/8c006b16dae02d30a174c0dc754488c542e54b64/data/cswa.rda), commit `8c006b16dae02d30a174c0dc754488c542e54b64`; [data documentation](https://github.com/cran/AHMbook/blob/8c006b16dae02d30a174c0dc754488c542e54b64/man/data_cswa.Rd). The original file, documentation, SHA-256 hashes and provenance are in `data/raw/`. Decoder used in this run: `rda_reader`.

- 43 points; 128 observed visits; 42 points with all 3 visits. Missing visits (`darty` visit 1) are excluded from the likelihood, not coded as 0.
- The source array has axes point × interval × visit; the code reshapes it to point × visit × interval.
- Plot areas are joined by point name. Points whose order differs between the count table and the covariate table: `zeala`, `zealb`, `zeal4a`.
- Interval lengths 2, 3, 5 minutes; first detections per interval [51, 14, 10], total 75 detections (not distinct individuals).
- Total plot area 23.37 ha; observed exposure 69.68 ha-visits.
- Spot-mapping: 38 patches, 74 individuals over 66.69 ha = **1.109612 individuals/ha** (descriptive only; not used in the likelihood).

## Primary analysis with an area offset

Set $\mu_i=\rho A_i$. For CRTE, the density column is the per-visit present density $\nu_i/A_i$.

| Model | logLik | AIC | ΔAIC | Weight | Density/ha | r per minute | c or φ |
| --- | --- | --- | --- | --- | --- | --- | --- |
| MPNc | -173.833004 | 353.666008 | 0.000000 | 0.950507 | 1.083596 | 0.373167 | 0.456632 |
| CRTE | -176.788168 | 359.576336 | 5.910329 | 0.049493 | 1.088156 | 0.373191 | 0.424957 |
| MPN0 | -189.690004 | 383.380008 | 29.714000 | 0.000000 | 1.086521 | 0.373191 | — |
| MPN1 | -203.614175 | 411.228350 | 57.562343 | 0.000000 | 2.031182 | 0.086074 | — |

MPNc estimates $\widehat\rho=1.083596$ individuals/ha, -2.34% relative to spot-mapping; $\widehat c=0.456632$; 10-minute cumulative detection probability 0.990635. CRTE: present density 1.088156, $\widehat\phi=0.424957$, superpopulation density 2.560623 (log CI 1.684–3.895); this is a different estimand. MPN1: density 2.031182 (+83.1% relative to spot-mapping), r = 0.086074.

## Confidence intervals

Primary CIs are built on the log (density) or logit (probability) scale and back-transformed, as described in the manuscript. `confidence_intervals.csv` also keeps natural-scale Wald intervals for comparison only.

| Model | Parameter | Estimate | Covariance | 95% CI (log/logit) |
| --- | --- | --- | --- | --- |
| MPN0 | rho | 1.087 | cluster | (0.773, 1.528) |
| MPN0 | r | 0.373 | cluster | (0.295, 0.458) |
| MPN1 | rho | 2.031 | model | (1.473, 2.801) |
| MPN1 | r | 0.086 | model | (0.062, 0.119) |
| MPNc | rho | 1.084 | model | (0.793, 1.481) |
| MPNc | r | 0.373 | model | (0.296, 0.457) |
| MPNc | c | 0.457 | model | (0.289, 0.635) |
| CRTE | rho | 1.088 | model | (0.801, 1.479) |
| CRTE | r | 0.373 | model | (0.296, 0.457) |
| CRTE | phi | 0.425 | model | (0.277, 0.587) |

MPN0 uses the point-clustered sandwich; the other models use the inverse observed information. For MPN0, the sandwich/information SE ratio is 1.504 for density and 1.016 for r: dependence between visits inflates mainly the variance of abundance, consistent with the Corollary on the rank-one correction (this is a descriptive observation; the design has an offset, so the balanced formula does not apply exactly).

## Removal composition test

Pooled over the observed visits, the interval counts are $\mathrm{Multinomial}(\sum d;\mathbf q(r))$ under every candidate model. Composition $\widehat r$ = 0.373191; expected (45.96, 22.42, 6.62) versus observed [51, 14, 10]. Pearson $X^2=5.444$ (p = 0.0196), $G^2=5.681$ (p = 0.0171), 1 degree of freedom; parametric bootstrap p for $G^2$ (4999 draws, seed 20261012) = 0.0222. The constant-hazard model fits the pooled composition poorly, suggesting heterogeneity in detectability.

## Third-order diagnostic

Area-based residual ratio over the 42 points with all three visits: 0.332202 using the MPNc estimates (predicted $\widehat c$ = 0.456632); 0.326303 using the CRTE estimates (predicted $\widehat\phi^2$ = 0.180589). Pooled version with per-visit means subtracted (audit only): 0.376415. This is a descriptive statistic, not a test.

## Sensitivity analysis: intercept only

The abundance column is in individuals/point/visit.

| Model | AIC | Weight | Individuals/point/visit | r per minute | c or φ |
| --- | --- | --- | --- | --- | --- |
| MPNc | 356.120210 | 0.836938 | 0.589969 | 0.371773 | 0.456554 |
| CRTE | 359.391472 | 0.163060 | 0.591729 | 0.373191 | 0.427994 |
| MPN0 | 381.743239 | 0.000002 | 0.591475 | 0.373191 | — |
| MPN1 | 410.844304 | 0.000000 | 1.100106 | 0.086601 | — |

## Numerical checks

50/50 optimizer runs converged; 8/8 fits have CIs. Largest log-likelihood spread across starting points 1.42e-13; largest |score| at the solution 1.72e-06; logLik change when tightening the tail truncation threshold from $10^{-12}$ to $10^{-14}$ and adding 16 terms: 0. These are numerical checks, not a proof of global optimality.

`results/paper_comparison.csv` lists every comparison with its rounding tolerance; `results/run_manifest.json` records software versions, configuration and hashes; `results/site_loglik.csv` stores per-point contributions.

## Scope

Points are treated as independent even though they belong to 38 patches. AIC weights are not probabilities that a biological mechanism is true; the data have no individual identities, so $c$ does not directly confirm that individuals return.
