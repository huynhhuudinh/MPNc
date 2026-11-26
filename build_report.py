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
                     ['Mô hình', 'logLik', 'AIC', 'ΔAIC', 'Trọng số', 'Mật độ/ha', 'r mỗi phút', 'c hoặc φ'])
    tab_sens = table([[f['model'], f"{f['aic']:.6f}", f"{f['akaike_weight']:.6f}", f"{f['estimate']:.6f}",
                       f"{f['r']:.6f}", dep(f)] for f in sens],
                     ['Mô hình', 'AIC', 'Trọng số', 'Số cá thể/điểm/lần', 'r mỗi phút', 'c hoặc φ'])
    prim = intervals[intervals.primary & (intervals.analysis == 'area')]
    tab_ci = table([[r.model, r.parameter, f"{r.estimate:.3f}", r.covariance,
                     f"({r.transformed_lower:.3f}, {r.transformed_upper:.3f})"]
                    for r in prim.itertuples()],
                   ['Mô hình', 'Tham số', 'Ước lượng', 'Hiệp phương sai', 'CI 95% (log/logit)'])
    m0 = fits['area', 'MPN0']
    ratio_rho = m0['mu_se_cluster'] / m0['mu_se_model']
    ratio_r = m0['r_se_cluster'] / m0['r_se_model']
    joins = ', '.join(f"`{m['counts_id']}`" for m in audit['alignment_mismatches_corrected'])
    c, crte, m1 = fits['area', 'MPNc'], fits['area', 'CRTE'], fits['area', 'MPN1']
    spot = audit['spot_density']
    status = ('Tất cả khớp.' if mismatches.empty else
              f"**{len(mismatches)} giá trị không khớp:** " +
              ', '.join(f"{r.analysis}/{r.model}/{r.quantity}" for r in mismatches.itertuples()))
    report = f'''# Kết quả tái lập phần Application

Đã chạy **8 phép khớp: 4 mô hình × 2 cách mô hình hóa abundance** trên dữ liệu thực. {diagnostics['fits_successful']}/{diagnostics['fits_total']} phép khớp thành công. Đã đối chiếu {len(comparison)} con số in trong Mục 7 của bản thảo (ước lượng, logLik, AIC, trọng số, khoảng tin cậy, chẩn đoán bậc ba, kiểm định composition, tổng diện tích): {int(comparison.matches_rounding.sum())}/{len(comparison)} khớp sau làm tròn. {status} {best} có AIC nhỏ nhất ở phân tích chính; {best_sens} ở phân tích chỉ có hệ số chặn.

## Dữ liệu

Nguồn: [AHMbook cswa.rda]({link}), commit `{source['commit']}`; [tài liệu dữ liệu]({doc}). Bản gốc, tài liệu, SHA-256 và nguồn nằm trong `data/raw/`. Bộ giải mã dùng ở lần chạy này: `{manifest.get('rda_decoder', 'rdata')}`.

- {audit['sites']} điểm; {audit['observed_visits']} lần khảo sát quan sát được; {audit['complete_sites']} điểm đủ 3 lần. Lần thiếu ({', '.join(f"`{m['point_id']}` lần {m['visit']}" for m in audit['missing_visits'])}) được loại khỏi likelihood, không mã hóa thành 0.
- Mảng nguồn có trục điểm × khoảng × lần; mã chuyển thành điểm × lần × khoảng.
- Diện tích được ghép theo tên điểm. Các điểm có thứ tự khác giữa bảng đếm và bảng covariate: {joins}.
- Khoảng dài {', '.join(f'{d:g}' for d in audit['intervals_minutes'])} phút; số phát hiện lần đầu theo khoảng {audit['counts_per_interval']}, tổng {audit['total_detections']} lượt phát hiện (không phải số cá thể riêng biệt).
- Tổng diện tích điểm {audit['sum_plot_area_ha']:.2f} ha; exposure quan sát {audit['observed_area_exposure_ha_visits']:.2f} ha-lần.
- Spot-mapping: {audit['spot_patches']} patch, {audit['spot_count']} cá thể trên {audit['spot_area_ha']:.2f} ha = **{spot:.6f} cá thể/ha** (chỉ dùng để mô tả, không vào likelihood).

## Phân tích chính với offset diện tích

Đặt $\\mu_i=\\rho A_i$. Với CRTE, cột mật độ là mật độ hiện diện trong một lần khảo sát $\\nu_i/A_i$.

{tab_main}

{c['model']} ước tính $\\widehat\\rho={c['estimate']:.6f}$ cá thể/ha, chênh {100*(c['estimate']/spot-1):+.2f}% so với spot-mapping; $\\widehat c={c['dep']:.6f}$; xác suất phát hiện tích lũy 10 phút {c['cumulative_detection']:.6f}. CRTE: mật độ hiện diện {crte['estimate']:.6f}, $\\widehat\\phi={crte['dep']:.6f}$, mật độ superpopulation {crte['superpopulation_estimate']:.6f} (CI log {crte['superpopulation_lo']:.3f}–{crte['superpopulation_hi']:.3f}); đây là estimand khác. MPN1: mật độ {m1['estimate']:.6f} ({100*(m1['estimate']/spot-1):+.1f}% so với spot-mapping), r = {m1['r']:.6f}.

## Khoảng tin cậy

CI chính được dựng trên thang log (mật độ) hoặc logit (xác suất) rồi biến đổi ngược, đúng như mô tả trong bản thảo. `confidence_intervals.csv` giữ thêm Wald trên thang gốc chỉ để đối chiếu.

{tab_ci}

MPN0 dùng sandwich gom cụm theo điểm; các mô hình khác dùng nghịch đảo observed information. Tỷ số SE sandwich/SE information của MPN0 là {ratio_rho:.3f} cho mật độ và {ratio_r:.3f} cho r: phụ thuộc giữa các lần làm tăng chủ yếu phương sai của abundance, phù hợp với Corollary về hiệu chỉnh hạng một (đây là quan sát mô tả; thiết kế có offset nên công thức cân bằng không áp dụng chính xác).

## Kiểm định thành phần removal

Gộp qua các lần quan sát, số đếm theo khoảng là $\\mathrm{{Multinomial}}(\\sum d;\\mathbf q(r))$ dưới mọi mô hình ứng viên. $\\widehat r$ từ composition = {gof['r_composition_mle']:.6f}; kỳ vọng ({', '.join(f'{e:.2f}' for e in gof['expected'])}) so với quan sát {gof['observed']}. Pearson $X^2={gof['pearson_x2']:.3f}$ (p = {gof['pearson_p']:.4f}), $G^2={gof['likelihood_ratio_g2']:.3f}$ (p = {gof['likelihood_ratio_p']:.4f}), {gof['df']} bậc tự do; p bootstrap tham số cho $G^2$ ({gof['bootstrap_draws']} lần, seed {gof['bootstrap_seed']}) = {gof['likelihood_ratio_bootstrap_p']:.4f}. Mô hình hazard hằng khớp kém với composition gộp, gợi ý dị biệt khả năng phát hiện.

## Chẩn đoán bậc ba

Tỷ số residual theo diện tích trên {third['complete_sites']} điểm đủ ba lần: {third['models']['MPNc']['area_residual_ratio']:.6f} khi dùng ước lượng MPNc (dự đoán $\\widehat c$ = {third['models']['MPNc']['predicted_ratio']:.6f}); {third['models']['CRTE']['area_residual_ratio']:.6f} khi dùng ước lượng CRTE (dự đoán $\\widehat\\phi^2$ = {third['models']['CRTE']['predicted_ratio']:.6f}). Bản gộp trừ trung bình theo lần (chỉ để audit): {third['models']['MPNc']['pooled_ratio_using_mean_mu']:.6f}. Đây là thống kê mô tả, không phải kiểm định.

## Phân tích độ nhạy: chỉ có hệ số chặn

Cột abundance có đơn vị cá thể/điểm/lần khảo sát.

{tab_sens}

## Kiểm tra số

{converged}/{attempts} lần tối ưu hội tụ; {diagnostics['intervals_available']}/{diagnostics['fits_total']} phép khớp có CI. Chênh log-likelihood lớn nhất giữa các điểm khởi đầu {maxspread:.3g}; |score| lớn nhất tại nghiệm {maxgrad:.3g}; thay đổi logLik khi siết ngưỡng cắt đuôi từ $10^{{-12}}$ xuống $10^{{-14}}$ và thêm 16 số hạng: {maxtail:.3g}. Đây là kiểm tra số, không phải chứng minh tối ưu toàn cục.

`results/paper_comparison.csv` liệt kê từng phép đối chiếu với dung sai làm tròn; `results/run_manifest.json` lưu phiên bản phần mềm, cấu hình và hash; `results/site_loglik.csv` lưu đóng góp theo điểm.

## Phạm vi

Các điểm được coi là độc lập dù thuộc {audit['spot_patches']} patch. Trọng số AIC không phải xác suất cơ chế sinh học đúng; dữ liệu không có định danh cá thể nên $c$ không xác nhận trực tiếp cá thể quay lại.
'''
    (ROOT / 'APPLICATION_REPORT.md').write_text(report, encoding='utf-8')
    print(json.dumps(diagnostics, indent=2))
    if not mismatches.empty:
        print(mismatches.to_string(index=False))


if __name__ == '__main__':
    main()
