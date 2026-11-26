#!/usr/bin/env python3
"""Build the application report and audit every Section 7 number of the manuscript.

All statements in the report are computed from the results folder; the only
hard-coded numbers are the manuscript's published values, which are used as
comparison targets and never as optimizer inputs.
"""
import argparse
import json
from pathlib import Path

import pandas as pd
from scipy.stats import norm

from run_application import ROOT

MODELS = ('MPN0', 'MPN1', 'MPNc', 'CRTE')

# Values printed in Section 7 of the revised manuscript (rounded to 3 decimals
# unless stated). Keys: (analysis, model, quantity).
PUBLISHED = {
    ('area', 'MPNc', 'loglik'): -173.833, ('area', 'MPNc', 'aic'): 353.666,
    ('area', 'MPNc', 'estimate'): 1.084, ('area', 'MPNc', 'r'): .373, ('area', 'MPNc', 'dep'): .457,
    ('area', 'CRTE', 'loglik'): -176.788, ('area', 'CRTE', 'aic'): 359.576,
    ('area', 'CRTE', 'estimate'): 1.088, ('area', 'CRTE', 'r'): .373, ('area', 'CRTE', 'dep'): .425,
    ('area', 'MPN0', 'loglik'): -189.690, ('area', 'MPN0', 'aic'): 383.380,
    ('area', 'MPN0', 'estimate'): 1.087, ('area', 'MPN0', 'r'): .373,
    ('area', 'MPN1', 'loglik'): -203.614, ('area', 'MPN1', 'aic'): 411.228,
    ('area', 'MPN1', 'estimate'): 2.031, ('area', 'MPN1', 'r'): .086,
    ('area', 'MPNc', 'delta_aic'): 0., ('area', 'CRTE', 'delta_aic'): 5.910,
    ('area', 'MPN0', 'delta_aic'): 29.714, ('area', 'MPN1', 'delta_aic'): 57.562,
    ('area', 'MPNc', 'akaike_weight'): .951, ('area', 'CRTE', 'akaike_weight'): .049,
    ('area', 'MPNc', 'cumulative_detection'): .990635, ('area', 'MPN1', 'cumulative_detection'): .593,
    ('area', 'CRTE', 'superpopulation_estimate'): 2.561,
    ('area', 'CRTE', 'superpopulation_lo'): 1.684, ('area', 'CRTE', 'superpopulation_hi'): 3.895,
    ('area', 'MPNc', 'mu_lo'): .793, ('area', 'MPNc', 'mu_hi'): 1.481,
    ('area', 'MPNc', 'r_lo'): .296, ('area', 'MPNc', 'r_hi'): .457,
    ('area', 'MPNc', 'dep_lo'): .289, ('area', 'MPNc', 'dep_hi'): .635,
    ('area', 'CRTE', 'mu_lo'): .801, ('area', 'CRTE', 'mu_hi'): 1.479,
    ('area', 'CRTE', 'dep_lo'): .277, ('area', 'CRTE', 'dep_hi'): .587,
    ('area', 'MPN0', 'mu_lo'): .773, ('area', 'MPN0', 'mu_hi'): 1.528,
    ('area', 'MPN0', 'r_lo'): .295, ('area', 'MPN0', 'r_hi'): .458,
    ('area', 'MPN1', 'mu_lo'): 1.473, ('area', 'MPN1', 'mu_hi'): 2.801,
    ('area', 'MPN1', 'r_lo'): .062, ('area', 'MPN1', 'r_hi'): .119,
    ('intercept', 'MPNc', 'aic'): 356.120, ('intercept', 'CRTE', 'aic'): 359.391,
    ('intercept', 'MPN0', 'aic'): 381.743, ('intercept', 'MPN1', 'aic'): 410.844,
    ('intercept', 'MPNc', 'akaike_weight'): .837, ('intercept', 'CRTE', 'akaike_weight'): .163,
    ('intercept', 'MPNc', 'estimate'): .590, ('intercept', 'CRTE', 'estimate'): .592,
    ('intercept', 'MPN0', 'estimate'): .591, ('intercept', 'MPN1', 'estimate'): 1.100,
    ('intercept', 'MPNc', 'r'): .372, ('intercept', 'CRTE', 'r'): .373,
    ('intercept', 'MPN0', 'r'): .373, ('intercept', 'MPN1', 'r'): .087,
    ('intercept', 'MPNc', 'dep'): .457, ('intercept', 'CRTE', 'dep'): .428,
    ('intercept', 'MPNc', 'mu_lo'): .432, ('intercept', 'MPNc', 'mu_hi'): .806,
    ('intercept', 'MPNc', 'dep_lo'): .287, ('intercept', 'MPNc', 'dep_hi'): .637,
}
PUBLISHED_OTHER = {
    ('third_order', 'MPNc', 'area_residual_ratio'): .332202,
    ('third_order', 'CRTE', 'area_residual_ratio'): .326303,
    ('composition', '-', 'r_composition_mle'): .373,
    ('composition', '-', 'pearson_x2'): 5.44, ('composition', '-', 'likelihood_ratio_g2'): 5.68,
    ('composition', '-', 'pearson_p'): .020, ('composition', '-', 'likelihood_ratio_p'): .017,
    ('data', '-', 'spot_density'): 1.109612, ('data', '-', 'sum_plot_area_ha'): 23.37,
    ('data', '-', 'observed_area_exposure_ha_visits'): 69.68,
}


def decimals(value):
    text = repr(float(value)).rstrip('0').split('.')
    return len(text[1]) if len(text) > 1 else 0


def table(rows, headers):
    return '\n'.join(['| ' + ' | '.join(headers) + ' |',
                      '| ' + ' | '.join(['---'] * len(headers)) + ' |'] +
                     ['| ' + ' | '.join(str(v) for v in row) + ' |' for row in rows])


def lookup(fits, analysis, model, key):
    f = fits[analysis, model]
    if key.endswith(('_lo', '_hi')) and not key.startswith('superpopulation'):
        return f[f"{key}_{f['primary_se']}"]
    return f[key]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--results', type=Path, default=ROOT / 'results')
    p = parser.parse_args().results
    manifest = json.loads((p / 'run_manifest.json').read_text())
    audit = manifest['data']
    fits = {(a, m): json.loads((p / f'{a}_{m}.json').read_text())
            for a in ('area', 'intercept') for m in MODELS}
    third = json.loads((p / 'third_order.json').read_text())
    gof = json.loads((p / 'composition_gof.json').read_text())

    # ---- audit against the manuscript --------------------------------------
    compare = []
    for (analysis, model, key), expected in PUBLISHED.items():
        actual = lookup(fits, analysis, model, key)
        compare.append((analysis, model, key, expected, actual))
    for (group, model, key), expected in PUBLISHED_OTHER.items():
        source = {'third_order': lambda: third['models'][model][key],
                  'composition': lambda: gof[key], 'data': lambda: audit[key]}[group]
        compare.append((group, model, key, expected, source()))
    rows = []
    for group, model, key, expected, actual in compare:
        tolerance = 0.5 * 10 ** -decimals(expected) + 1e-9
        rows.append(dict(analysis=group, model=model, quantity=key, paper=expected,
                         reproduced=actual, difference=actual - expected,
                         rounding_tolerance=tolerance,
                         matches_rounding=abs(actual - expected) <= tolerance))
    comparison = pd.DataFrame(rows)
    comparison.to_csv(p / 'paper_comparison.csv', index=False)
    mismatches = comparison[~comparison.matches_rounding]

    # ---- interval table --------------------------------------------------------
    interval_rows = []
    z = norm.ppf(.975)
    for (analysis, model), f in fits.items():
        params = [('mu', f['estimate']), ('r', f['r'])]
        if model in ('MPNc', 'CRTE'):
            params.append(('dep', f['dep']))
        for param, estimate in params:
            for method in ('model', 'cluster'):
                se = f[f'{param}_se_{method}']
                name = (('rho' if analysis == 'area' else 'mu') if param == 'mu' else
                        ('c' if model == 'MPNc' else 'phi') if param == 'dep' else 'r')
                interval_rows.append(dict(
                    analysis=analysis, model=model, parameter=name, covariance=method,
                    primary=(method == f['primary_se']), estimate=estimate, se_natural=se,
                    transformed_lower=f[f'{param}_lo_{method}'],
                    transformed_upper=f[f'{param}_hi_{method}'],
                    natural_wald_lower=estimate - z * se, natural_wald_upper=estimate + z * se))
    intervals = pd.DataFrame(interval_rows)
    intervals.to_csv(p / 'confidence_intervals.csv', index=False)

    # ---- diagnostics -----------------------------------------------------------
    maxgrad = max(f.get('score_max_abs', 0.) for f in fits.values())
    maxspread = max(f['starts_nll_spread'] for f in fits.values())
    maxtail = max(f['truncation_loglik_change'] for f in fits.values())
    converged = sum(t['optimizer_success'] for f in fits.values() for t in f['attempts'])
    attempts = sum(len(f['attempts']) for f in fits.values())
    diagnostics = dict(
        fits_successful=sum(f['success'] for f in fits.values()), fits_total=len(fits),
        intervals_available=sum(f['interval_valid'] for f in fits.values()),
        optimizer_attempts=attempts, optimizer_attempts_converged=converged,
        maximum_absolute_score=maxgrad, maximum_start_nll_spread=maxspread,
        maximum_truncation_loglik_change=maxtail,
        manuscript_comparisons=len(comparison),
        manuscript_comparisons_matching_rounding=int(comparison.matches_rounding.sum()))
    (p / 'verification_summary.json').write_text(json.dumps(diagnostics, indent=2) + '\n')

    # ---- report ----------------------------------------------------------------
    source = manifest['provenance']
    doc = f"https://github.com/cran/AHMbook/blob/{source['commit']}/man/data_cswa.Rd"
    link = f"https://github.com/cran/AHMbook/blob/{source['commit']}/data/cswa.rda"
    primary = sorted([f for (a, m), f in fits.items() if a == 'area'], key=lambda f: f['aic'])
    sens = sorted([f for (a, m), f in fits.items() if a == 'intercept'], key=lambda f: f['aic'])
    best, best_sens = primary[0]['model'], sens[0]['model']
    dep = lambda f: f"{f['dep']:.6f}" if f['model'] in ('MPNc', 'CRTE') else '—'
    tab_main = table([[f['model'], f"{f['loglik']:.6f}", f"{f['aic']:.6f}", f"{f['delta_aic']:.6f}",
                       f"{f['akaike_weight']:.6f}", f"{f['estimate']:.6f}", f"{f['r']:.6f}", dep(f)]
                      for f in primary],
                     ['Model', 'logLik', 'AIC', 'ΔAIC', 'Weight', 'Density/ha', 'r per minute', 'c or φ'])
    tab_sens = table([[f['model'], f"{f['aic']:.6f}", f"{f['akaike_weight']:.6f}", f"{f['estimate']:.6f}",
                       f"{f['r']:.6f}", dep(f)] for f in sens],
                     ['Model', 'AIC', 'Weight', 'Individuals/point/visit', 'r per minute', 'c or φ'])
    prim = intervals[intervals.primary & (intervals.analysis == 'area')]
    tab_ci = table([[r.model, r.parameter, f"{r.estimate:.3f}", r.covariance,
                     f"({r.transformed_lower:.3f}, {r.transformed_upper:.3f})"]
                    for r in prim.itertuples()],
                   ['Model', 'Parameter', 'Estimate', 'Covariance', '95% CI (log/logit)'])
    m0 = fits['area', 'MPN0']
    ratio_rho = m0['mu_se_cluster'] / m0['mu_se_model']
    ratio_r = m0['r_se_cluster'] / m0['r_se_model']
    joins = ', '.join(f"`{m['counts_id']}`" for m in audit['alignment_mismatches_corrected'])
    c, crte, m1 = fits['area', 'MPNc'], fits['area', 'CRTE'], fits['area', 'MPN1']
    spot = audit['spot_density']
    status = ('All match.' if mismatches.empty else
              f"**{len(mismatches)} values do not match:** " +
              ', '.join(f"{r.analysis}/{r.model}/{r.quantity}" for r in mismatches.itertuples()))
    report = f'''# Application reproduction results

Ran **8 fits: 4 models × 2 abundance specifications** on the real data. {diagnostics['fits_successful']}/{diagnostics['fits_total']} fits succeeded. Checked {len(comparison)} numbers printed in Section 7 of the manuscript (estimates, logLik, AIC, weights, confidence intervals, third-order diagnostic, composition test, total area): {int(comparison.matches_rounding.sum())}/{len(comparison)} match after rounding. {status} {best} has the smallest AIC in the primary analysis; {best_sens} in the intercept-only analysis.

## Data

Source: [AHMbook cswa.rda]({link}), commit `{source['commit']}`; [data documentation]({doc}). The original file, documentation, SHA-256 hashes and provenance are in `data/raw/`. Decoder used in this run: `{manifest.get('rda_decoder', 'rdata')}`.

- {audit['sites']} points; {audit['observed_visits']} observed visits; {audit['complete_sites']} points with all 3 visits. Missing visits ({', '.join(f"`{m['point_id']}` visit {m['visit']}" for m in audit['missing_visits'])}) are excluded from the likelihood, not coded as 0.
- The source array has axes point × interval × visit; the code reshapes it to point × visit × interval.
- Plot areas are joined by point name. Points whose order differs between the count table and the covariate table: {joins}.
- Interval lengths {', '.join(f'{d:g}' for d in audit['intervals_minutes'])} minutes; first detections per interval {audit['counts_per_interval']}, total {audit['total_detections']} detections (not distinct individuals).
- Total plot area {audit['sum_plot_area_ha']:.2f} ha; observed exposure {audit['observed_area_exposure_ha_visits']:.2f} ha-visits.
- Spot-mapping: {audit['spot_patches']} patches, {audit['spot_count']} individuals over {audit['spot_area_ha']:.2f} ha = **{spot:.6f} individuals/ha** (descriptive only; not used in the likelihood).

## Primary analysis with an area offset

Set $\\mu_i=\\rho A_i$. For CRTE, the density column is the per-visit present density $\\nu_i/A_i$.

{tab_main}

{c['model']} estimates $\\widehat\\rho={c['estimate']:.6f}$ individuals/ha, {100*(c['estimate']/spot-1):+.2f}% relative to spot-mapping; $\\widehat c={c['dep']:.6f}$; 10-minute cumulative detection probability {c['cumulative_detection']:.6f}. CRTE: present density {crte['estimate']:.6f}, $\\widehat\\phi={crte['dep']:.6f}$, superpopulation density {crte['superpopulation_estimate']:.6f} (log CI {crte['superpopulation_lo']:.3f}–{crte['superpopulation_hi']:.3f}); this is a different estimand. MPN1: density {m1['estimate']:.6f} ({100*(m1['estimate']/spot-1):+.1f}% relative to spot-mapping), r = {m1['r']:.6f}.

## Confidence intervals

Primary CIs are built on the log (density) or logit (probability) scale and back-transformed, as described in the manuscript. `confidence_intervals.csv` also keeps natural-scale Wald intervals for comparison only.

{tab_ci}

MPN0 uses the point-clustered sandwich; the other models use the inverse observed information. For MPN0, the sandwich/information SE ratio is {ratio_rho:.3f} for density and {ratio_r:.3f} for r: dependence between visits inflates mainly the variance of abundance, consistent with the Corollary on the rank-one correction (this is a descriptive observation; the design has an offset, so the balanced formula does not apply exactly).

## Removal composition test

Pooled over the observed visits, the interval counts are $\\mathrm{{Multinomial}}(\\sum d;\\mathbf q(r))$ under every candidate model. Composition $\\widehat r$ = {gof['r_composition_mle']:.6f}; expected ({', '.join(f'{e:.2f}' for e in gof['expected'])}) versus observed {gof['observed']}. Pearson $X^2={gof['pearson_x2']:.3f}$ (p = {gof['pearson_p']:.4f}), $G^2={gof['likelihood_ratio_g2']:.3f}$ (p = {gof['likelihood_ratio_p']:.4f}), {gof['df']} degree{'' if gof['df'] == 1 else 's'} of freedom; parametric bootstrap p for $G^2$ ({gof['bootstrap_draws']} draws, seed {gof['bootstrap_seed']}) = {gof['likelihood_ratio_bootstrap_p']:.4f}. The constant-hazard model fits the pooled composition poorly, suggesting heterogeneity in detectability.

## Third-order diagnostic

Area-based residual ratio over the {third['complete_sites']} points with all three visits: {third['models']['MPNc']['area_residual_ratio']:.6f} using the MPNc estimates (predicted $\\widehat c$ = {third['models']['MPNc']['predicted_ratio']:.6f}); {third['models']['CRTE']['area_residual_ratio']:.6f} using the CRTE estimates (predicted $\\widehat\\phi^2$ = {third['models']['CRTE']['predicted_ratio']:.6f}). Pooled version with per-visit means subtracted (audit only): {third['models']['MPNc']['pooled_ratio_using_mean_mu']:.6f}. This is a descriptive statistic, not a test.

## Sensitivity analysis: intercept only

The abundance column is in individuals/point/visit.

{tab_sens}

## Numerical checks

{converged}/{attempts} optimizer runs converged; {diagnostics['intervals_available']}/{diagnostics['fits_total']} fits have CIs. Largest log-likelihood spread across starting points {maxspread:.3g}; largest |score| at the solution {maxgrad:.3g}; logLik change when tightening the tail truncation threshold from $10^{{-12}}$ to $10^{{-14}}$ and adding 16 terms: {maxtail:.3g}. These are numerical checks, not a proof of global optimality.

`results/paper_comparison.csv` lists every comparison with its rounding tolerance; `results/run_manifest.json` records software versions, configuration and hashes; `results/site_loglik.csv` stores per-point contributions.

## Scope

Points are treated as independent even though they belong to {audit['spot_patches']} patches. AIC weights are not probabilities that a biological mechanism is true; the data have no individual identities, so $c$ does not directly confirm that individuals return.
'''
    (ROOT / 'APPLICATION_REPORT.md').write_text(report, encoding='utf-8')
    print(json.dumps(diagnostics, indent=2))
    if not mismatches.empty:
        print(mismatches.to_string(index=False))


if __name__ == '__main__':
    main()
