#!/usr/bin/env python3
"""Fit the paper's four removal models to the public AHMbook cswa data.

Raw R data are preserved and hash-checked. No manuscript estimate initializes
the optimizer. Missing visits are integrated out; areas are already hectares.
"""
from __future__ import annotations

from softclosure import (FitConfig, NumericalFitError, Likelihood,
                         binomial_table, intervals, numerical_derivatives,
                         stable_gradient, write_json)
import argparse
from dataclasses import asdict
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import platform
import time
import urllib.request

import numpy as np
import pandas as pd
import scipy
from scipy.optimize import minimize, minimize_scalar
from scipy.special import expit, logit, gammaln, logsumexp, xlogy
from scipy.stats import poisson, norm

ROOT = Path(__file__).resolve().parent
MODELS = ('MPN0', 'MPN1', 'MPNc', 'CRTE')
DURATIONS = np.array([2., 3., 5.])


def unequal_probabilities(r, durations=DURATIONS):
    durations = np.asarray(durations, float)
    if not 0 < r < 1 or durations.ndim != 1 or np.any(durations <= 0):
        raise ValueError('Invalid detection probability or interval durations')
    log_survival = np.log1p(-r)
    starts = np.r_[0., np.cumsum(durations)[:-1]]
    pi = np.exp(starts * log_survival) * -np.expm1(durations * log_survival)
    a = float(-np.expm1(durations.sum() * log_survival))
    return pi, a, pi / a


def verify_sources(data_dir, download=False):
    raw = Path(data_dir) / 'raw'
    provenance = json.loads((raw / 'provenance.json').read_text())
    for item in provenance['files']:
        target = raw / item['filename']
        if download:
            with urllib.request.urlopen(item['url'], timeout=60) as response:
                content = response.read()
            if sha256(content).hexdigest() != item['sha256']:
                raise ValueError(f"Download hash mismatch: {target.name}")
            target.write_bytes(content)
        if not target.exists() or sha256(target.read_bytes()).hexdigest() != item['sha256']:
            raise ValueError(f"Missing or changed source file: {target}; use --download")
    return provenance


def read_source(path):
    """Decode cswa.rda; points x intervals x occasions, as stored by R.

    Uses the ``rdata`` package when installed and otherwise the dependency-free
    ``rda_reader`` module. Both return the same arrays (checked in the tests).
    """
    try:
        import rdata
    except ImportError:
        from rda_reader import read_cswa
        values, dimnames, covs, spot = read_cswa(path)
        return values, dimnames, covs, spot, 'rda_reader'
    obj = rdata.read_rda(path)['cswa']
    counts = obj['counts']
    dimnames = [[str(v) for v in counts.coords[d].values] for d in counts.dims]
    return (np.asarray(counts.values, float), dimnames, obj['covs'].copy(),
            obj['spotMaps'].copy(), f'rdata {rdata.__version__}')


def load_data(data_dir=ROOT / 'data', export=False):
    """R dimensions: points x intervals x occasions, NOT points x occasions x intervals."""
    data_dir = Path(data_dir)
    verify_sources(data_dir)
    values, dimnames, covs, spot, _ = read_source(data_dir / 'raw' / 'cswa.rda')
    ids = dimnames[0]
    if len(set(ids)) != len(ids) or not covs.index.is_unique or set(ids) != set(covs.index):
        raise ValueError('Nonunique or inconsistent point identifiers')
    if dimnames[1] != ['int1', 'int2', 'int3']:
        raise ValueError('Unexpected interval axis')
    if dimnames[2] != ['occ1', 'occ2', 'occ3']:
        raise ValueError('Unexpected occasion axis')
    positional_mismatches = [dict(position=i+1, counts_id=a, covariates_id=str(b))
                             for i, (a, b) in enumerate(zip(ids, covs.index)) if a != b]
    covs = covs.loc[ids]  # Explicit identifier join, never positional assignment.
    y = np.asarray(values, float).transpose(0, 2, 1)
    area = covs['plotArea'].to_numpy(dtype=float)  # Source documentation: ha.
    like = ApplicationLikelihood(y, area)
    missing = [dict(point_id=ids[i], visit=int(t+1)) for i, t in zip(*np.where(~like.observed))]
    audit = dict(sites=len(ids), visits_per_site=3, intervals_minutes=DURATIONS.tolist(),
                 observed_visits=int(like.observed.sum()), complete_sites=int(like.observed.all(axis=1).sum()),
                 missing_visits=missing, count_axis_original=['point','interval','occasion'],
                 count_axis_analysis=['point','occasion','interval'],
                 alignment_mismatches_corrected=positional_mismatches,
                 counts_per_interval=np.nansum(y, axis=(0,1)).astype(int).tolist(),
                 total_detections=int(np.nansum(y)), maximum_visit_total=like.maxd,
                 plot_area_units='ha', sum_plot_area_ha=float(area.sum()),
                 observed_area_exposure_ha_visits=float(np.dot(area, like.nvisits)),
                 spot_patches=len(spot), spot_count=int(spot.CSWA.sum()),
                 spot_area_ha=float(spot.parea.sum()),
                 spot_density=float(spot.CSWA.sum()/spot.parea.sum()))
    if export:
        records=[]
        for i, point in enumerate(ids):
            for t in range(y.shape[1]):
                for j, duration in enumerate(DURATIONS):
                    records.append(dict(point_id=point, visit=t+1, interval=j+1,
                                        duration_minutes=duration, count=y[i,t,j],
                                        observed=bool(like.observed[i,t]), plot_area_ha=area[i]))
        pd.DataFrame(records).to_csv(data_dir/'counts_long.csv', index=False)
        covs.rename_axis('point_id').to_csv(data_dir/'site_covariates.csv')
        spot.rename_axis('patch_id').to_csv(data_dir/'spot_mapping.csv')
        np.savez_compressed(data_dir/'analysis_data.npz', y=y, area_ha=area, point_id=ids)
        write_json(data_dir/'data_audit.json', audit)
    return y, area, ids, audit


class ApplicationLikelihood(Likelihood):
    """Area-dependent site means, unequal intervals, and fully missing visits.

    A common integer support covers every site's Poisson prior. Guard terms
    protect relative log accuracy in addition to the absolute prior-tail cutoff.
    """
    def __init__(self, y, exposure, config=None, durations=DURATIONS):
        y = np.asarray(y, dtype=float)
        if y.ndim != 3 or not y.size:
            raise ValueError('Require nonempty sites x visits x intervals array')
        partial = np.isnan(y).any(axis=2) != np.isnan(y).all(axis=2)
        if partial.any():
            raise ValueError('Partially missing interval counts require a different observation likelihood')
        observed_values = y[~np.isnan(y)]
        if np.any(~np.isfinite(observed_values)) or np.any(observed_values < 0) or np.any(observed_values != np.floor(observed_values)):
            raise ValueError('Observed counts must be finite nonnegative integers')
        self.y = y.copy()
        self.S, self.T, self.J = y.shape
        self.observed = ~np.isnan(y).all(axis=2)
        self.nvisits = self.observed.sum(axis=1)
        if np.any(self.nvisits == 0):
            raise ValueError('Every included site must have at least one observed visit')
        self.exposure = np.asarray(exposure, dtype=float)
        if self.exposure.shape != (self.S,) or np.any(~np.isfinite(self.exposure)) or np.any(self.exposure <= 0):
            raise ValueError('Each site requires a finite positive exposure')
        self.durations = np.asarray(durations, dtype=float)
        if self.durations.shape != (self.J,) or np.any(self.durations <= 0):
            raise ValueError('Duration vector must match intervals')
        self.cfg = config or FitConfig()
        # Zero below is only an indexing placeholder; the observed mask removes
        # every missing visit's likelihood contribution in all four models.
        filled = np.nan_to_num(y, nan=0.).astype(int)
        self.d = filled.sum(axis=2)
        self.maxd = int(self.d.max())
        self.counts = filled.sum(axis=1)
        self.logfactorials = gammaln(filled + 1).sum(axis=(1,2))
        self.composition_constants = gammaln(self.d+1).sum(axis=1)-self.logfactorials

    def site_loglik(self, model, mu, r, dep=None, tail=None, pad=0):
        if not np.isfinite(mu) or mu <= 0:
            raise NumericalFitError('Invalid abundance/density')
        pi, a, q = unequal_probabilities(r, self.durations)
        means = mu * self.exposure
        composition = self.composition_constants + self.counts @ np.log(q)
        if model == 'MPN0' or (model in ('MPNc','CRTE') and dep == 0):
            return (self.counts.sum(axis=1)*np.log(means) + self.counts@np.log(pi)
                    - self.nvisits*means*a - self.logfactorials)
        if model in ('MPNc','CRTE') and dep == 1:
            model = 'MPN1'
        if model == 'MPN1':
            lam, p = means, a
        elif model in ('MPNc','CRTE') and dep is not None and 0 < dep < 1:
            lam, p = (means*dep, None) if model == 'MPNc' else (means/dep, dep*a)
        else:
            raise NumericalFitError('Invalid model or dependence parameter')
        threshold = self.cfg.tail if tail is None else tail
        guard = max(8, self.maxd) + pad
        low = max(0, int(poisson.ppf(threshold/2, lam.min()))-guard)
        high = max(self.maxd, int(poisson.isf(threshold/2, lam.max()))+guard)
        if high-low+1 > self.cfg.support_cap:
            raise NumericalFitError('Latent support exceeds configured memory cap')
        k = np.arange(low, high+1)
        lp = -lam[:,None] + xlogy(k[None,:], lam[:,None])-gammaln(k[None,:]+1)
        if model == 'MPNc':
            b = binomial_table(k, self.maxd, a)
            residual = (1-dep)*means*a
            zz = np.arange(self.maxd+1)
            poisson_part = (-residual[:,None] + xlogy(zz[None,:],residual[:,None])
                            - gammaln(zz[None,:]+1))
            conditional = np.empty((self.S,self.maxd+1,len(k)))
            for d in range(self.maxd+1):
                conditional[:,d,:] = logsumexp(
                    b[None,:d+1,:] + poisson_part[:,d::-1,None], axis=1)
            for t in range(self.T):
                term = conditional[np.arange(self.S), self.d[:,t], :]
                lp += np.where(self.observed[:,t,None], term, 0.)
        else:
            conditional = binomial_table(k, self.maxd, p)
            for t in range(self.T):
                lp += np.where(self.observed[:,t,None], conditional[self.d[:,t]], 0.)
        result = composition + logsumexp(lp, axis=1)
        if not np.isfinite(result).all():
            raise NumericalFitError('Nonfinite site log likelihood')
        return result


def working_start(like):
    counts = like.counts.sum(axis=0)
    if counts.sum() == 0:
        raise NumericalFitError('All counts are zero: no positive interior abundance MLE')
    res = minimize_scalar(lambda eta: -float(counts@np.log(
        unequal_probabilities(float(expit(eta)),like.durations)[2])),
        bounds=(logit(like.cfg.r_min),logit(like.cfg.r_max)), method='bounded',
        options={'xatol':1e-12})
    r = float(expit(res.x))
    mu = float(counts.sum()/(np.dot(like.exposure,like.nvisits)*unequal_probabilities(r,like.durations)[1]))
    return mu, r


def fit_one(like, model, starts=8, prior_fits=None):
    begin = time.perf_counter()
    mu0, r0 = working_start(like)
    cfg = like.cfg
    bounds=[(np.log(cfg.mu_min),np.log(cfg.mu_max)),(logit(cfg.r_min),logit(cfg.r_max))]
    if model in ('MPNc','CRTE'):
        bounds.append((logit(cfg.dep_min),logit(cfg.dep_max)))
    guesses=[]
    # Deterministic dispersed starts using only data, not manuscript estimates.
    settings=[(r0,.2),(r0,.5),(r0,.8),(.08,.3),(.18,.7),(.60,.5),(.25,.9),(.45,.1)]
    if model == 'MPN0':
        settings=[(r0,None)]
    for r, dep in settings[:starts]:
        mu = mu0*unequal_probabilities(r0,like.durations)[1]/unequal_probabilities(r,like.durations)[1]
        guesses.append(np.array([np.log(mu),logit(r)]+([logit(dep)] if len(bounds)==3 else [])))
    while len(guesses)<starts and model != 'MPN0':
        rng=np.random.default_rng(700+len(guesses))
        guesses.append(guesses[0]+rng.normal(0,.6,len(bounds)))
    fun=lambda x:like.objective(model,x)
    attempts=[]
    for initial in guesses:
        initial=np.clip(initial,np.array(bounds)[:,0]+1e-6,np.array(bounds)[:,1]-1e-6)
        opt=minimize(fun,initial,jac=lambda x:stable_gradient(fun,x),method='L-BFGS-B',bounds=bounds,
                     options={'maxiter':cfg.maxiter,'ftol':2e-13,'gtol':1e-6,'maxls':40})
        gradient=float(np.max(np.abs(stable_gradient(fun,opt.x))))
        if not opt.success or gradient>1e-3:
            pre=minimize(fun,opt.x,method='Powell',bounds=bounds,
                         options={'maxiter':cfg.maxiter,'ftol':1e-12,'xtol':1e-8})
            polished=minimize(fun,pre.x if pre.fun<opt.fun else opt.x,
                              jac=lambda x:stable_gradient(fun,x),method='L-BFGS-B',bounds=bounds,
                              options={'maxiter':cfg.maxiter*2,'ftol':1e-14,'gtol':1e-6,'maxls':50})
            if polished.fun<=opt.fun+1e-9:
                opt=polished
        attempts.append(dict(initial_x=initial.tolist(),x=opt.x.tolist(),nll=float(opt.fun),
                             optimizer_success=bool(opt.success),message=str(opt.message),
                             gradient_max_abs=float(np.max(np.abs(stable_gradient(fun,opt.x))))))
    valid=[v for v in attempts if np.isfinite(v['nll']) and v['nll']<1e90]
    if not valid:
        raise NumericalFitError(f'No finite candidate for {model}')
    best=min(valid,key=lambda v:v['nll'])
    x=np.array(best['x'])
    mu,r=float(np.exp(x[0])),float(expit(x[1]))
    dep=float(expit(x[2])) if len(x)==3 else (0. if model=='MPN0' else 1.)
    boundary=None
    boundary_checks=[]
    if len(x)==3:
        for other,d in [('MPN0',0.),('MPN1',1.)]:
            f=(prior_fits or {}).get(other)
            if f and f['success']:
                boundary_checks.append(dict(model=other,loglik=f['loglik']))
                if f['loglik']>=-best['nll']-1e-8:
                    mu,r,dep=f['estimate'],f['r'],d
                    boundary=('independence_limit' if model=='CRTE' and d==0 else f'dep={d:g}')
                    x=np.array([np.log(mu),logit(r)])
    base=like.site_loglik(model,mu,r,dep)
    tight=like.site_loglik(model,mu,r,dep,tail=cfg.check_tail,pad=16)
    bound_hit=any(min(v-lo,hi-v)<2e-5 for v,(lo,hi) in zip(best['x'],bounds)) if boundary is None else False
    result=dict(model=model,estimate=mu,r=r,dep=dep,x=x.tolist(),boundary=boundary,
                loglik=float(base.sum()),k=len(bounds),aic=2*len(bounds)-2*float(base.sum()),
                cumulative_detection=unequal_probabilities(r,like.durations)[1],
                truncation_loglik_change=float(abs(tight.sum()-base.sum())),
                truncation_max_site_change=float(np.max(np.abs(tight-base))),
                attempts=attempts,boundary_checks=boundary_checks,
                starts_nll_spread=float(max(v['nll'] for v in valid)-min(v['nll'] for v in valid)),
                optimization_bound=bool(bound_hit),interval_valid=False,
                primary_se='cluster' if model=='MPN0' else 'model')
    result['success']=bool(best['optimizer_success'] and not bound_hit and
                           (boundary is not None or best['gradient_max_abs']<5e-3) and
                           result['truncation_loglik_change']<cfg.likelihood_tolerance)
    if boundary is None and result['success']:
        ci=intervals(like,model,x)
        result.update(ci)
        result['interval_valid']=True
        # Save both covariance matrices for derived estimands and reproducibility.
        H,sc=numerical_derivatives(like,model,x,cfg.hessian_step)
        cov=np.linalg.inv(H)
        result['covariance_x_model']=cov.tolist()
        result['covariance_x_cluster']=(cov@(sc.T@sc)@cov).tolist()
        if model=='CRTE':
            # log(lambda_density) = log(nu_density) - log(phi)
            g=np.array([1.,0.,-(1-dep)])
            se=float(np.sqrt(g@cov@g))
            lam=mu/dep
            result.update(superpopulation_estimate=lam,superpopulation_lo=float(lam*np.exp(-norm.ppf(.975)*se)),
                          superpopulation_hi=float(lam*np.exp(norm.ppf(.975)*se)))
    result['seconds']=time.perf_counter()-begin
    return result


def _composition_r(counts, durations, cfg):
    res = minimize_scalar(lambda eta: -float(counts @ np.log(
        unequal_probabilities(float(expit(eta)), durations)[2])),
        bounds=(logit(cfg.r_min), logit(cfg.r_max)), method='bounded',
        options={'xatol': 1e-12})
    return float(expit(res.x))


def composition_gof(like, bootstrap=4999, seed=20261012):
    """Goodness of fit of the constant-hazard removal composition.

    Conditional on the detected totals, every candidate model gives independent
    Multinomial(d_it; q(r)) compositions with a common q(r). The pooled interval
    counts are therefore Multinomial(sum d; q(r)) whatever the latent dependence,
    and this check does not depend on the abundance or dependence construction.
    """
    from scipy.stats import chi2
    counts = like.counts.sum(axis=0).astype(float)
    total = float(counts.sum())
    _, r = working_start(like)  # maximizes counts @ log q(r)
    expected = total * unequal_probabilities(r, like.durations)[2]
    df = like.J - 2
    if df < 1:
        raise ValueError('Composition goodness of fit needs at least three intervals')
    pearson = float(np.sum((counts - expected) ** 2 / expected))
    deviance = float(2 * np.sum(xlogy(counts, counts / expected)))
    # Parametric bootstrap of G^2 under the fitted composition, refitting r each
    # time, checks the chi-square approximation with small expected counts.
    rng = np.random.default_rng(seed)
    q_hat = unequal_probabilities(r, like.durations)[2]
    exceed = 0
    for draw in rng.multinomial(int(total), q_hat, size=bootstrap):
        draw = draw.astype(float)
        e = total * unequal_probabilities(_composition_r(draw, like.durations, like.cfg),
                                          like.durations)[2]
        exceed += 2 * np.sum(xlogy(draw, draw / e)) >= deviance - 1e-12
    boot_p = float((exceed + 1) / (bootstrap + 1))
    return dict(observed=counts.astype(int).tolist(), expected=expected.tolist(),
                bootstrap_draws=bootstrap, bootstrap_seed=seed, likelihood_ratio_bootstrap_p=boot_p,
                total_detections=int(total), r_composition_mle=r, df=df,
                pearson_x2=pearson, pearson_p=float(chi2.sf(pearson, df)),
                likelihood_ratio_g2=deviance, likelihood_ratio_p=float(chi2.sf(deviance, df)),
                note='Pooled over observed visits; valid under every candidate because '
                     'compositions are conditionally independent of latent dependence.')


def third_order(y, area, fits):
    """Descriptive diagnostics; no hypothesis test or bootstrap claim.

    Primary ratio uses area-dependent fitted means and is model-specific because
    rho and r differ by fit. Pooled centering is supplied only for audit.
    """
    complete=~np.isnan(y).any(axis=(1,2))
    d=y[complete].sum(axis=2)
    n=len(d)
    pooled=float(np.mean(np.prod(d-d.mean(axis=0),axis=1)))
    out={'complete_sites':n,'pooled_third_central_moment':pooled,'models':{}}
    for model in ('MPNc','CRTE'):
        f=fits[model];a=f['cumulative_detection'];mu=f['estimate']*area[complete]
        residual_product=np.prod(d-mu[:,None]*a,axis=1)
        ratio=float(residual_product.sum()/np.sum(mu*a**3))
        out['models'][model]=dict(area_residual_ratio=ratio,
            predicted_ratio=f['dep'] if model=='MPNc' else f['dep']**2,
            pooled_ratio_using_mean_mu=float(pooled/(mu.mean()*a**3)))
    return out


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir',type=Path,default=ROOT/'data')
    parser.add_argument('--out',type=Path,default=ROOT/'results')
    parser.add_argument('--download',action='store_true',help='Redownload the pinned source files and verify SHA-256')
    parser.add_argument('--prepare-only',action='store_true')
    parser.add_argument('--starts',type=int,default=8)
    args=parser.parse_args()
    if args.starts<1:
        parser.error('--starts must be positive')
    provenance=verify_sources(args.data_dir,args.download)
    y,area,ids,audit=load_data(args.data_dir,export=True)
    print(json.dumps(audit,indent=2),flush=True)
    if args.prepare_only:
        return
    args.out.mkdir(parents=True,exist_ok=True)
    cfg=FitConfig(starts=args.starts,maxiter=350)
    all_results={};rows=[];site_rows=[]
    for analysis,exposure in [('area',area),('intercept',np.ones(len(area)))]:
        like=ApplicationLikelihood(y,exposure,cfg)
        fits={}
        for model in MODELS:
            print(f'Fitting {analysis}/{model} ...',flush=True)
            f=fit_one(like,model,args.starts,fits)
            f['analysis']=analysis
            f['estimand']='visit_density_birds_per_ha' if analysis=='area' else 'visit_abundance_birds_per_point'
            fits[model]=f
            write_json(args.out/f'{analysis}_{model}.json',f)
            print(f"  logLik={f['loglik']:.9f}; AIC={f['aic']:.9f}; estimate={f['estimate']:.9f}; r={f['r']:.9f}; dep={f['dep']:.9f}; success={f['success']}",flush=True)
        aics=np.array([fits[m]['aic'] for m in MODELS])
        delta=aics-aics.min();weights=np.exp(-delta/2);weights/=weights.sum()
        for model,dd,ww in zip(MODELS,delta,weights):
            f=fits[model];f.update(delta_aic=float(dd),akaike_weight=float(ww))
            write_json(args.out/f'{analysis}_{model}.json',f)
            row={k:v for k,v in f.items() if not isinstance(v,(list,dict))}
            rows.append(row)
            sl=like.site_loglik(model,f['estimate'],f['r'],f['dep'])
            for point,ll in zip(ids,sl):
                site_rows.append(dict(analysis=analysis,model=model,point_id=point,loglik=ll))
        all_results[analysis]=fits
    pd.DataFrame(rows).to_csv(args.out/'model_results.csv',index=False)
    pd.DataFrame(site_rows).to_csv(args.out/'site_loglik.csv',index=False)
    write_json(args.out/'third_order.json',third_order(y,area,all_results['area']))
    write_json(args.out/'composition_gof.json',composition_gof(ApplicationLikelihood(y,area,cfg)))
    write_json(args.out/'run_manifest.json',dict(
        completed_utc=datetime.now(timezone.utc).isoformat(),python=platform.python_version(),
        numpy=np.__version__,scipy=scipy.__version__,pandas=pd.__version__,
        rda_decoder=read_source(args.data_dir/'raw'/'cswa.rda')[4],
        config=asdict(cfg),data=audit,provenance=provenance,
        source_hashes={p.name:sha256(p.read_bytes()).hexdigest() for p in ROOT.glob('*.py')},
        all_fits_successful=all(f['success'] for group in all_results.values() for f in group.values())))
    if not all(f['success'] for group in all_results.values() for f in group.values()):
        raise SystemExit('One or more fits failed; inspect diagnostic JSON files.')


if __name__=='__main__':
    main()
