#!/usr/bin/env python3
"""Reproducible repeated-removal simulations; NumPy/SciPy only for inference.

Models: MPN0, MPN1, MPNc, CRTE. See README.md for estimands and limitations.
All fits include the full multinomial composition, not merely visit totals.
"""
from __future__ import annotations

import os
for _key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
             "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_key] = "1"

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from hashlib import sha256
from itertools import product, combinations
import json
from pathlib import Path
import platform
import sys
import time
import warnings

import numpy as np
import scipy
from scipy.optimize import brentq, minimize
from scipy.special import expit, logit, gammaln, logsumexp, xlogy
from scipy.stats import poisson, norm

MODELS = ("MPN0", "MPN1", "MPNc", "CRTE")
VERSION = "1.0.0"


@dataclass(frozen=True)
class Scenario:
    block: str
    dgm: str
    S: int
    T: int
    mu: float
    r: float
    dep: float
    J: int = 3

    @property
    def id(self):
        def f(x):
            return f"{x:g}".replace(".", "p")
        return (f"{self.block}_{self.dgm}_S{self.S}_T{self.T}_J{self.J}"
                f"_mu{f(self.mu)}_r{f(self.r)}_c{f(self.dep)}")


@dataclass(frozen=True)
class FitConfig:
    starts: int = 3
    maxiter: int = 250
    tail: float = 1e-12
    check_tail: float = 1e-14
    support_cap: int = 12000
    mu_min: float = 1e-4
    mu_max: float = 500.0
    r_min: float = 1e-5
    r_max: float = 1 - 1e-5
    dep_min: float = 1e-4
    dep_max: float = 1 - 1e-4
    hessian_step: float = 2e-4
    intervals: bool = True
    likelihood_tolerance: float = 1e-6
    tie_tolerance: float = 1e-6


class NumericalFitError(RuntimeError):
    pass


def pass_probabilities(r, J):
    pi = r * np.power(1.0 - r, np.arange(J))
    a = -np.expm1(J * np.log1p(-r))
    return pi, float(a), pi / a


def seeded_rng(base, scenario_id, replicate, stream=0):
    # Stable across Python versions, process scheduling and subsets of scenarios.
    words = np.frombuffer(sha256(scenario_id.encode()).digest()[:16], dtype="<u4")
    return np.random.default_rng(np.random.SeedSequence(
        [int(base), int(replicate), int(stream), *map(int, words)]))


def generate(scenario: Scenario, rng):
    s = scenario
    if s.dgm not in ("MPNc", "CRTE") or not (s.mu > 0 and 0 < s.r < 1
                                             and 0 <= s.dep <= 1):
        raise ValueError("Invalid generating parameters")
    if s.dgm == "MPNc":
        common = rng.poisson(s.dep * s.mu, size=(s.S, 1))
        N = common + rng.poisson((1 - s.dep) * s.mu, size=(s.S, s.T))
    else:
        if s.dep <= 0:
            raise ValueError("CRTE generation requires positive availability")
        M = rng.poisson(s.mu / s.dep, size=(s.S, 1))
        N = rng.binomial(np.broadcast_to(M, (s.S, s.T)), s.dep)
    remaining = N.copy()
    y = np.empty((s.S, s.T, s.J), dtype=np.int64)
    # Sequential binomial removals exactly reproduce the multinomial kernel.
    for j in range(s.J):
        y[:, :, j] = rng.binomial(remaining, s.r)
        remaining -= y[:, :, j]
    return y


def poisson_support(lam, tail, cap, pad=0, minimum_upper=0):
    if not np.isfinite(lam) or lam < 0:
        raise NumericalFitError("Invalid mixing mean")
    if lam == 0:
        return np.array([0], dtype=int), np.array([0.0]), 0.0
    lo = max(0, int(poisson.ppf(tail / 2, lam)) - pad)
    hi = max(minimum_upper, int(poisson.isf(tail / 2, lam)) + pad)
    if hi - lo + 1 > cap:
        raise NumericalFitError(f"Latent support width exceeds {cap}")
    k = np.arange(lo, hi + 1)
    logp = -lam + xlogy(k, lam) - gammaln(k + 1)
    omitted = float(poisson.cdf(lo - 1, lam) + poisson.sf(hi, lam))
    return k, logp, omitted


def binomial_table(k, maxd, p):
    """log P(Bin(k,p)=d), rows d and columns k; no renormalization."""
    d = np.arange(maxd + 1)[:, None]
    kk = np.asarray(k)[None, :]
    valid = kk >= d
    rem = np.maximum(kk - d, 0)
    ans = (gammaln(kk + 1) - gammaln(d + 1) - gammaln(rem + 1)
           + xlogy(d, p) + xlogy(rem, 1 - p))
    return np.where(valid, ans, -np.inf)


def conditional_total_table(k, maxd, mu, r, c, J):
    """Exact Bin(k,a)+Pois((1-c)*mu*a) convolution, log scale."""
    _, a, _ = pass_probabilities(r, J)
    b = binomial_table(k, maxd, a)
    w = (1 - c) * mu * a
    z = np.arange(maxd + 1)
    lp = -w + xlogy(z, w) - gammaln(z + 1)
    out = np.empty_like(b)
    for d in range(maxd + 1):
        out[d] = logsumexp(b[:d + 1] + lp[d::-1, None], axis=0)
    return out


class Likelihood:
    def __init__(self, y, config=None):
        y = np.asarray(y)
        if y.ndim != 3 or np.any(y < 0) or np.any(y != np.floor(y)):
            raise ValueError("y must be a nonnegative integer S x T x J array")
        self.y = y.astype(np.int64)
        self.S, self.T, self.J = y.shape
        self.cfg = config or FitConfig()
        self.d = self.y.sum(axis=2)
        self.patterns, self.inverse, self.frequency = np.unique(
            self.d, axis=0, return_inverse=True, return_counts=True)
        self.maxd = int(self.d.max())
        self.counts = self.y.sum(axis=1)
        self.composition_constants = (gammaln(self.d + 1).sum(axis=1)
                                      - gammaln(self.y + 1).sum(axis=(1, 2)))

    def site_loglik(self, model, mu, r, dep=None, tail=None, pad=0):
        if mu <= 0 or not 0 < r < 1:
            raise NumericalFitError("Parameters outside natural domain")
        pi, a, q = pass_probabilities(r, self.J)
        composition = self.composition_constants + self.counts @ np.log(q)
        if model == "MPN0" or (model in ("MPNc", "CRTE") and dep == 0):
            return (-mu * pi[None, None, :] + self.y * np.log(mu * pi)
                    - gammaln(self.y + 1)).sum(axis=(1, 2))
        if model in ("MPNc", "CRTE") and dep == 1:
            model = "MPN1"
        threshold = self.cfg.tail if tail is None else tail
        if model == "MPN1":
            lam, p = mu, a
        elif model == "CRTE":
            if dep is None or not 0 < dep < 1:
                raise NumericalFitError("Invalid CRTE availability")
            lam, p = mu / dep, dep * a
        elif model == "MPNc":
            if dep is None or not 0 < dep < 1:
                raise NumericalFitError("Invalid soft-closure parameter")
            lam, p = dep * mu, None
        else:
            raise ValueError(model)
        # Prior tail alone need not control relative probability at an unlikely
        # observation. Add an observation-aware guard on BOTH support ends.
        guard = max(8, self.maxd)
        k, logprior, _ = poisson_support(lam, threshold, self.cfg.support_cap,
                                        pad + guard, self.maxd if p is not None else 0)
        if model == "MPNc":
            conditional = conditional_total_table(k, self.maxd, mu, r, dep, self.J)
        else:
            conditional = binomial_table(k, self.maxd, p)
        # Matrix has one row per distinct total-count pattern, one column per k.
        joint = logprior[None, :] + conditional[self.patterns].sum(axis=1)
        lp = logsumexp(joint, axis=1)
        result = composition + lp[self.inverse]
        if not np.all(np.isfinite(result)):
            raise NumericalFitError("Nonfinite site log likelihood")
        return result

    def from_x(self, model, x, tail=None, pad=0):
        mu, r = float(np.exp(x[0])), float(expit(x[1]))
        dep = float(expit(x[2])) if model in ("MPNc", "CRTE") else None
        return self.site_loglik(model, mu, r, dep, tail, pad)

    def objective(self, model, x):
        try:
            return -float(self.from_x(model, x).sum())
        except (NumericalFitError, ValueError, FloatingPointError):
            return 1e100


def working_fit(like):
    cfg = like.cfg
    totals = like.y.sum(axis=(0, 1))
    M = float(totals.sum())
    if M == 0:
        raise NumericalFitError("All-zero sample: positive abundance MLE absent")
    mean_h = float(np.dot(np.arange(like.J), totals) / M)
    if like.J < 2:
        raise NumericalFitError("Working fit needs at least two removal passes")

    def eh(r):
        return float(np.dot(np.arange(like.J), pass_probabilities(r, like.J)[2]))

    if mean_h >= eh(cfg.r_min):
        r, edge = cfg.r_min, True
    elif mean_h <= eh(cfg.r_max):
        r, edge = cfg.r_max, True
    else:
        r, edge = brentq(lambda r: eh(r) - mean_h, cfg.r_min, cfg.r_max,
                         xtol=1e-13), False
    mu = M / (like.S * like.T * pass_probabilities(r, like.J)[1])
    edge = bool(edge or not cfg.mu_min < mu < cfg.mu_max)
    return dict(mu=float(mu), r=float(r), bound_hit=edge,
                x=np.array([np.log(mu), logit(r)]))


def initial_dependence(like, mu, r):
    if like.T < 2:
        return 0.5
    covariance = np.cov(like.d, rowvar=False, ddof=1)
    mean_cov = np.mean(covariance[np.triu_indices(like.T, 1)])
    a = pass_probabilities(r, like.J)[1]
    return float(np.clip(mean_cov / (mu * a * a), 0.02, 0.98))


def numerical_derivatives(like, model, x, step):
    """Hessian of total negative log likelihood and site log-score matrix."""
    x = np.asarray(x, float)
    n = len(x)
    h = step * np.maximum(1.0, np.abs(x))
    center = like.from_x(model, x, tail=like.cfg.check_tail)
    f0 = -center.sum()
    H = np.zeros((n, n))
    scores = np.empty((like.S, n))
    for j in range(n):
        e = np.zeros(n); e[j] = h[j]
        plus = like.from_x(model, x + e, tail=like.cfg.check_tail)
        minus = like.from_x(model, x - e, tail=like.cfg.check_tail)
        scores[:, j] = (plus - minus) / (2 * h[j])
        H[j, j] = (-plus.sum() - minus.sum() - 2 * f0) / h[j] ** 2
        for k in range(j):
            e2 = np.zeros(n); e2[k] = h[k]
            vals = [-like.from_x(model, x + a * e + b * e2,
                                 tail=like.cfg.check_tail).sum()
                    for a, b in ((1, 1), (1, -1), (-1, 1), (-1, -1))]
            H[j, k] = H[k, j] = (vals[0] - vals[1] - vals[2] + vals[3]) / (
                4 * h[j] * h[k])
    return H, scores


def intervals(like, model, x):
    H, scores = numerical_derivatives(like, model, x, like.cfg.hessian_step)
    eig = np.linalg.eigvalsh(H)
    if not np.all(np.isfinite(eig)) or eig.min() <= 1e-7:
        raise NumericalFitError("Observed information is not positive definite")
    condition = float(eig.max() / eig.min())
    if condition > 1e10:
        raise NumericalFitError("Ill-conditioned observed information")
    cov_model = np.linalg.inv(H)
    cov_cluster = cov_model @ (scores.T @ scores) @ cov_model
    z = norm.ppf(0.975)
    estimates = np.r_[np.exp(x[0]), expit(x[1:])]
    jac = np.r_[estimates[0], estimates[1:] * (1 - estimates[1:])]
    result = {"information_condition": condition,
              "score_max_abs": float(np.max(np.abs(scores.sum(axis=0))))}
    for name, cov in (("model", cov_model), ("cluster", cov_cluster)):
        if np.any(np.diag(cov) <= 0) or not np.all(np.isfinite(cov)):
            raise NumericalFitError("Invalid covariance")
        se = np.sqrt(np.diag(cov))
        lower_x, upper_x = np.asarray(x) - z * se, np.asarray(x) + z * se
        lower = np.r_[np.exp(lower_x[0]), expit(lower_x[1:])]
        upper = np.r_[np.exp(upper_x[0]), expit(upper_x[1:])]
        if not np.all(np.isfinite(lower)) or not np.all(np.isfinite(upper)):
            raise NumericalFitError("Interval endpoints overflowed")
        for j, param in enumerate(("mu", "r", "dep")[:len(x)]):
            result[f"{param}_se_{name}"] = float(jac[j] * se[j])
            result[f"{param}_lo_{name}"] = float(lower[j])
            result[f"{param}_hi_{name}"] = float(upper[j])
    return result


def blank_fit(model, message):
    return {"model": model, "success": False, "interval_valid": False,
            "error": str(message), "mu": None, "r": None, "dep": None,
            "loglik": None, "aic": None, "boundary": None}


def stable_gradient(fun, x):
    x = np.asarray(x, float)
    h = 1e-5 * np.maximum(1., np.abs(x))
    result = np.empty(len(x))
    for j in range(len(x)):
        e = np.zeros(len(x)); e[j] = h[j]
        result[j] = (fun(x + e) - fun(x - e)) / (2*h[j])
    return result


def fit_model(like, model, working, closed=None, rng=None):
    start_time = time.perf_counter()
    cfg = like.cfg
    mu0, r0 = working["mu"], working["r"]
    r_start = 0.15 if working["bound_hit"] else r0
    mu_start = float(like.d.mean()) / pass_probabilities(r_start, like.J)[1]
    c0 = initial_dependence(like, mu_start, r_start)
    bounds = [(np.log(cfg.mu_min), np.log(cfg.mu_max)),
              (logit(cfg.r_min), logit(cfg.r_max))]
    attempts = []
    boundary = None
    if model == "MPN0":
        x = working["x"].copy()
        mu, r, dep = mu0, r0, 0.0
        success = not working["bound_hit"]
        objective = like.objective(model, x)
        source = "analytic_working_profile"
        bound_hit = working["bound_hit"]
    else:
        if model == "MPN1":
            a0 = pass_probabilities(r_start, like.J)[1]
            effective_r = 1 - (1 - c0 * a0) ** (1 / like.J)
            starts = [[np.log(mu_start), logit(r_start)],
                      [np.log(min(mu_start / c0, cfg.mu_max / 2)), logit(effective_r)]]
        else:
            bounds += [(logit(cfg.dep_min), logit(cfg.dep_max))]
            starts = [[np.log(mu_start), logit(r_start), logit(c)]
                      for c in (c0, 0.25, 0.75, 0.5)]
        rng = rng or np.random.default_rng(123)
        while len(starts) < cfg.starts:
            starts.append((np.asarray(starts[0]) + rng.normal(0, 0.35, len(bounds))).tolist())
        for initial in starts[:cfg.starts]:
            initial = np.clip(initial, np.array(bounds)[:, 0] + 1e-7,
                              np.array(bounds)[:, 1] - 1e-7)
            fun = lambda v: like.objective(model, v)
            res = minimize(fun, initial, jac=lambda v: stable_gradient(fun, v),
                           method="L-BFGS-B", bounds=bounds,
                           options={"maxiter": cfg.maxiter, "ftol": 2e-12,
                                    "gtol": 2e-5, "maxls": 30})
            attempts.append(res)
        admissible = [v for v in attempts if np.isfinite(v.fun) and v.fun < 1e90]
        if not admissible:
            return blank_fit(model, "No finite optimization result")
        best = min(admissible, key=lambda v: v.fun)
        if not best.success or np.max(np.abs(stable_gradient(fun, best.x))) > 1e-3:
            # A small function change can trigger false convergence on a flat
            # likelihood. Polish that candidate with a different search method.
            pre = minimize(fun, best.x, method="Powell", bounds=bounds,
                           options={"maxiter": cfg.maxiter, "ftol": 1e-11, "xtol": 1e-7})
            retry = minimize(fun, pre.x if pre.fun < best.fun else best.x,
                             jac=lambda v: stable_gradient(fun, v), method="L-BFGS-B",
                             bounds=bounds, options={"maxiter": cfg.maxiter*2,
                             "ftol": 2e-13, "gtol": 1e-5, "maxls": 40})
            attempts.append(retry)
            if np.isfinite(retry.fun) and retry.fun < 1e90:
                admissible.append(retry)
                best = min(admissible, key=lambda v: v.fun)
        x = best.x.copy()
        objective = float(best.fun)
        source = "interior"
        mu, r = float(np.exp(x[0])), float(expit(x[1]))
        dep = 1.0 if model == "MPN1" else float(expit(x[2]))
        bound_hit = bool(any(min(v - lo, hi - v) < 2e-5
                             for v, (lo, hi) in zip(x, bounds)))
        stationarity = float(np.max(np.abs(stable_gradient(fun, x))))
        success = bool(best.success and not bound_hit and stationarity <= 5e-3)
        if model in ("MPNc", "CRTE"):
            # Evaluate exact boundary/limiting models independently; never mimic
            # c=0 or c=1 by extreme logits and call them ordinary interior fits.
            low_obj = -float(like.site_loglik(model, mu0, r0, 0).sum())
            candidates = [(objective, "interior", mu, r, dep, success)]
            candidates.append((low_obj, "independence_limit" if model == "CRTE"
                                else "c=0", mu0, r0, 0.0, not working["bound_hit"]))
            if closed is not None and closed.get("success"):
                candidates.append((-closed["loglik"], "c=1" if model == "MPNc"
                                   else "phi=1", closed["mu"], closed["r"],
                                   1.0, True))
            # Select boundary in an objective tie; record that it is nonregular.
            minimum = min(v[0] for v in candidates)
            candidates = [v for v in candidates if v[0] <= minimum + 1e-8]
            pick = min(candidates, key=lambda v: v[1] == "interior")
            objective, source, mu, r, dep, selected_success = pick
            if source != "interior":
                boundary = source
                x = np.array([np.log(mu), logit(r)])
                bound_hit = False
                # A failed interior search is retained as a failure to establish
                # the global candidate comparison, even if a boundary is finite.
                success = bool(selected_success and any(v.success for v in admissible))
            else:
                success = selected_success
    result = dict(model=model, mu=mu, r=r, dep=dep, success=bool(success),
                  boundary=boundary, optimization_bound=bound_hit, source=source,
                  loglik=-float(objective), aic=2 * (2 if model in ("MPN0", "MPN1") else 3)
                  + 2 * float(objective), k=2 if model in ("MPN0", "MPN1") else 3,
                  interval_valid=False, error=None, primary_se="cluster" if model == "MPN0"
                  else "model", near_boundary=bool(dep is not None and model in ("MPNc", "CRTE")
                                                    and (dep < 0.02 or dep > 0.98)),
                  attempts=len(attempts), converged_attempts=sum(v.success for v in attempts),
                  best_start_nll_spread=float(max(v.fun for v in attempts if v.fun < 1e90)
                                             - min(v.fun for v in attempts if v.fun < 1e90))
                  if attempts and any(v.fun < 1e90 for v in attempts) else 0.0)
    try:
        base = like.site_loglik(model, mu, r, dep)
        tight = like.site_loglik(model, mu, r, dep, tail=cfg.check_tail, pad=16)
        result["truncation_loglik_change"] = float(abs(tight.sum() - base.sum()))
        result["truncation_max_site_change"] = float(np.max(np.abs(tight - base)))
        if result["truncation_loglik_change"] > cfg.likelihood_tolerance:
            result.update(success=False, error="Latent truncation sensitivity exceeded tolerance")
        if cfg.intervals and result["success"] and boundary is None:
            try:
                result.update(intervals(like, model, x))
                result["interval_valid"] = True
                if result["score_max_abs"] > 5e-3:
                    result.update(success=False, interval_valid=False,
                                  error="Score stationarity criterion exceeded tolerance")
            except (NumericalFitError, np.linalg.LinAlgError, ValueError) as exc:
                result["interval_error"] = str(exc)
        elif boundary is not None:
            result["interval_error"] = "Full-model boundary optimum: no ordinary Wald interval"
        if not result["success"] and result["error"] is None:
            result["error"] = ("Working profile has no admissible interior maximum" if model == "MPN0"
                               else "Optimizer, stationarity, or computational-bound criterion failed")
    except (NumericalFitError, ValueError, OverflowError) as exc:
        result.update(success=False, error=str(exc))
    result["seconds"] = time.perf_counter() - start_time
    return result


def fit_all(y, config=None, rng=None):
    like = Likelihood(y, config)
    try:
        working = working_fit(like)
    except NumericalFitError as exc:
        return [blank_fit(m, exc) for m in MODELS]
    result = []
    for model in MODELS:
        try:
            closed = result[1] if len(result) > 1 else None
            result.append(fit_model(like, model, working, closed, rng))
        except (NumericalFitError, ValueError, np.linalg.LinAlgError, OverflowError) as exc:
            result.append(blank_fit(model, exc))
    return result


def third_order(y, scenario, fits):
    d = y.sum(axis=2).astype(float)
    S, T = d.shape
    if T < 3 or S < 3:
        return {}
    centered = d - d.mean(axis=0)
    third = [np.mean(centered[:, a] * centered[:, b] * centered[:, c])
             for a, b, c in combinations(range(T), 3)]
    # Unbiased third cross-cumulant at each distinct triplet.
    k3 = float(np.mean(third) * S * S / ((S - 1) * (S - 2)))
    a = pass_probabilities(scenario.r, scenario.J)[1]
    k3true = scenario.mu * a ** 3 * (scenario.dep if scenario.dgm == "MPNc"
                                    else scenario.dep ** 2)
    out = {"k3_empirical": k3, "k3_true": k3true}
    ref = fits[0]
    if ref.get("success"):
        denominator = ref["mu"] * pass_probabilities(ref["r"], scenario.J)[1] ** 3
        out["r3_empirical"] = k3 / denominator
        out["r3_true_reference"] = scenario.dep if scenario.dgm == "MPNc" else scenario.dep ** 2
        for fit in fits[2:]:
            if fit.get("success"):
                a_fit = pass_probabilities(fit["r"], scenario.J)[1]
                pred = fit["mu"] * a_fit ** 3 * (fit["dep"] if fit["model"] == "MPNc"
                                                else fit["dep"] ** 2)
                out[f"k3_predicted_{fit['model']}"] = pred
                out[f"r3_predicted_{fit['model']}"] = pred / denominator
    return out


def scenarios(suite):
    core = [Scenario("core", "MPNc", 200, T, mu, r, c)
            for T, mu, r, c in product((3, 7), (2., 5.), (.25, .5), (.2, .5, .8))]
    reverse = [Scenario("reverse", "CRTE", v.S, v.T, v.mu, v.r, v.dep)
               for v in core]
    small = [Scenario("small", "MPNc", 43, 3, mu, .25, c)
             for mu, c in product((2., 5.), (.2, .5, .8))]
    stress = [Scenario("stress", "MPNc", 200, 3, 2., r, c)
              for r, c in ((.1, .2), (.1, .8), (.25, .05), (.25, .95), (.25, 0.), (.25, 1.))]
    if suite == "smoke":
        return [core[0], core[12], core[-1], small[1], reverse[0], reverse[-1],
                stress[0], stress[-2], stress[-1]]
    return {"core": core, "reverse": reverse, "small": small, "stress": stress,
            "full": core + reverse + small + stress}[suite]


def code_digest():
    return sha256(Path(__file__).read_bytes()).hexdigest()


def sanitize(obj):
    if isinstance(obj, dict):
        return {k: sanitize(v) for k, v in obj.items()}
    if isinstance(obj, (tuple, list)):
        return [sanitize(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return sanitize(obj.tolist())
    if isinstance(obj, (np.floating, float)):
        return float(obj) if np.isfinite(obj) else None
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, np.bool_):
        return bool(obj)
    return obj


def write_json(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(sanitize(obj), ensure_ascii=False, indent=2, allow_nan=False))
    temp.replace(path)


def task_fingerprint(scenario, cfg, base, digest):
    return sha256(json.dumps([asdict(scenario), asdict(cfg), base, digest],
                             sort_keys=True).encode()).hexdigest()


def run_one(scenario, replicate, base, cfg, digest):
    y = generate(scenario, seeded_rng(base, scenario.id, replicate))
    fits = fit_all(y, cfg, seeded_rng(base, scenario.id, replicate, stream=1))
    return sanitize({"scenario": asdict(scenario), "scenario_id": scenario.id,
                     "replicate": replicate, "fingerprint": task_fingerprint(scenario, cfg, base, digest),
                     "data_sha256": sha256(y.tobytes()).hexdigest(), "fits": fits,
                     "diagnostic": third_order(y, scenario, fits)})


def run_suite(suite, out, reps, workers, base, cfg):
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    selected = scenarios(suite)
    digest = code_digest()
    repdir = out / "replicates"; repdir.mkdir(exist_ok=True)
    manifest = {"version": VERSION, "code_sha256": digest, "suite": suite,
                "replicates_per_scenario": reps, "workers": workers, "base_seed": base,
                "fit_config": asdict(cfg), "scenarios": [dict(id=s.id, **asdict(s)) for s in selected],
                "created_utc": datetime.now(timezone.utc).isoformat(),
                "python": sys.version, "platform": platform.platform(),
                "numpy": np.__version__, "scipy": scipy.__version__,
                "purpose": "smoke_test_not_full_evidence" if suite == "smoke" else "simulation_run"}
    jobs = []
    retained = []
    for s in selected:
        expected = task_fingerprint(s, cfg, base, digest)
        for rep in range(reps):
            path = repdir / f"{s.id}_rep{rep:04d}.json"
            if path.exists():
                item = json.loads(path.read_text())
                if item.get("fingerprint") != expected:
                    raise RuntimeError(f"Resume fingerprint mismatch: {path}. Use a new output directory.")
                retained.append(item)
            else:
                jobs.append((s, rep, path))
    write_json(out / "manifest.json", manifest)
    print(f"{suite}: {len(selected)} scenarios x {reps} replicates; {len(jobs)} tasks pending; "
          f"{len(retained)} resumed; {workers} workers", flush=True)
    start = time.perf_counter()
    with ProcessPoolExecutor(max_workers=workers) as pool:
        pending = {pool.submit(run_one, s, rep, base, cfg, digest): path for s, rep, path in jobs}
        for done, future in enumerate(as_completed(pending), 1):
            record = future.result()  # Programming errors abort; they are not silently censored.
            write_json(pending[future], record)
            retained.append(record)
            if done == 1 or done % max(1, len(jobs) // 12) == 0 or done == len(jobs):
                print(f"Completed {done}/{len(jobs)}; elapsed {time.perf_counter()-start:.1f}s", flush=True)
    manifest["elapsed_seconds_this_invocation"] = time.perf_counter() - start
    manifest["completed_datasets"] = len(retained)
    manifest["completed_fits"] = sum(len(x["fits"]) for x in retained)
    write_json(out / "manifest.json", manifest)
    summarize(out, retained)
    return retained


def safe_mcse(values):
    v = np.asarray(values, float)
    return float(v.std(ddof=1) / np.sqrt(len(v))) if len(v) > 1 else np.nan


def summarize(out, records=None):
    import pandas as pd
    out = Path(out)
    if records is None:
        manifest = json.loads((out / "manifest.json").read_text())
        records = []
        for raw_scenario in manifest["scenarios"]:
            sid = raw_scenario["id"]
            scenario = Scenario(**{k: v for k, v in raw_scenario.items() if k != "id"})
            expected = task_fingerprint(scenario, FitConfig(**manifest["fit_config"]),
                                        manifest["base_seed"], manifest["code_sha256"])
            for rep in range(manifest["replicates_per_scenario"]):
                path = out / "replicates" / f"{sid}_rep{rep:04d}.json"
                if path.exists():
                    item = json.loads(path.read_text())
                    if item.get("fingerprint") != expected:
                        raise RuntimeError(f"Report fingerprint mismatch: {path}")
                    records.append(item)
        if not records:
            raise RuntimeError("No completed records for this manifest")
    rows, diagnostic_rows, selection_rows = [], [], []
    for item in records:
        s = item["scenario"]
        metadata = {"scenario_id": item["scenario_id"], "replicate": item["replicate"],
                    "block": s["block"], "dgm": s["dgm"], "S": s["S"], "T": s["T"],
                    "J": s["J"], "true_mu": s["mu"], "true_r": s["r"], "true_dep": s["dep"]}
        for fit in item["fits"]:
            rows.append(dict(metadata, **fit))
        diagnostic_rows.append(dict(metadata, **item["diagnostic"]))
        bymodel = {f["model"]: f for f in item["fits"]}
        m, c = bymodel["MPNc"], bymodel["CRTE"]
        pair_ok = m["success"] and c["success"]
        delta = c["aic"] - m["aic"] if pair_ok else np.nan
        identifiable = s["T"] >= 3 and 0 < s["dep"] < 1
        tie = bool(pair_ok and abs(delta) <= 1e-6)
        correct = np.nan
        if pair_ok and identifiable:
            correct = 0.5 if tie else float((delta > 0) == (s["dgm"] == "MPNc"))
        all_ok = all(f["success"] for f in bymodel.values())
        sel = dict(metadata, pair_valid=pair_ok, delta_aic_crte_minus_mpnc=delta,
                   pair_tie=tie, mechanisms_distinct=identifiable, correct_selection_credit=correct,
                   decisive_correct=float((delta > 2) if s["dgm"] == "MPNc" else (delta < -2))
                   if pair_ok and identifiable else np.nan, all_four_valid=all_ok)
        for name in MODELS:
            sel[f"winner_credit_{name}"] = np.nan
        if all_ok:
            minimum = min(f["aic"] for f in bymodel.values())
            winners = [name for name, f in bymodel.items() if f["aic"] <= minimum + 1e-6]
            for name in MODELS:
                sel[f"winner_credit_{name}"] = 1 / len(winners) if name in winners else 0.
        selection_rows.append(sel)
    fits = pd.DataFrame(rows).sort_values(["scenario_id", "replicate", "model"])
    fits.to_csv(out / "fits.csv", index=False)
    summaries = []
    groups = ["scenario_id", "block", "dgm", "S", "T", "J", "true_mu", "true_r", "true_dep", "model"]
    for key, g in fits.groupby(groups, dropna=False, sort=True):
        meta = dict(zip(groups, key)); model = meta["model"]
        for param in (("mu", "r", "dep") if model in ("MPNc", "CRTE") else ("mu", "r")):
            true = meta[f"true_{param}"]
            valid = g[g.success & g[param].notna()]
            values = valid[param].to_numpy(float)
            errors = values - true
            n = len(values); N = len(g)
            rmse = np.sqrt(np.mean(errors ** 2)) if n else np.nan
            primary = "cluster" if model == "MPN0" else "model"
            lower, upper = f"{param}_lo_{primary}", f"{param}_hi_{primary}"
            ci = valid.dropna(subset=[lower, upper]) if lower in valid and upper in valid else valid.iloc[:0]
            if param == "dep" and true in (0., 1.):
                ci = ci.iloc[:0]  # Do not call boundary Wald coverage nominal coverage.
            covered = ((ci[lower] <= true) & (true <= ci[upper])).to_numpy(float) if len(ci) else np.array([])
            status = ("generating_parameter" if param != "dep" or model == meta["dgm"]
                      else "second_order_reference")
            row = dict(meta, parameter=param, target=true, target_kind=status,
                       attempts=N, successful_fits=n, fit_success_rate=n / N,
                       convergence_mcse=np.sqrt((n/N)*(1-n/N)/N),
                       mean=np.mean(values) if n else np.nan,
                       median=np.median(values) if n else np.nan,
                       bias=np.mean(errors) if n else np.nan,
                       relative_bias_percent=100*np.mean(errors)/true if n and true != 0 else np.nan,
                       bias_mcse=safe_mcse(errors), empirical_sd=np.std(values, ddof=1) if n > 1 else np.nan,
                       scaled_mad=1.4826*np.median(np.abs(values-np.median(values))) if n else np.nan,
                       rmse=rmse, rmse_mcse=safe_mcse(errors**2)/(2*rmse) if n > 1 and rmse > 0 else np.nan,
                       interval_method=primary, valid_intervals=len(ci),
                       interval_availability=len(ci)/N,
                       coverage_given_valid=np.mean(covered) if len(ci) else np.nan,
                       coverage_mcse=np.sqrt(np.mean(covered)*(1-np.mean(covered))/len(ci)) if len(ci) else np.nan,
                       covered_per_attempt=np.sum(covered)/N if not
                       (param == "dep" and true in (0., 1.)) else np.nan,
                       mean_se=ci[f"{param}_se_{primary}"].mean() if len(ci) else np.nan,
                       mean_interval_width=(ci[upper]-ci[lower]).mean() if len(ci) else np.nan,
                       exact_boundary_rate=g.boundary.notna().mean(),
                       median_seconds=g.seconds.median() if "seconds" in g else np.nan)
            summaries.append(row)
    summary = pd.DataFrame(summaries)
    summary.to_csv(out / "parameter_summary.csv", index=False)
    selection = pd.DataFrame(selection_rows)
    selection.to_csv(out / "model_selection.csv", index=False)
    srows = []
    for sid, g in selection.groupby("scenario_id", sort=True):
        usable = g[g.pair_valid]
        credit = usable.correct_selection_credit.dropna().to_numpy()
        row = dict(scenario_id=sid, dgm=g.dgm.iloc[0], attempts=len(g),
                   valid_pairs=len(usable), paired_fit_availability=len(usable)/len(g),
                   correct_selection_rate=np.mean(credit) if len(credit) else np.nan,
                   selection_mcse=safe_mcse(credit),
                   median_delta_aic=usable.delta_aic_crte_minus_mpnc.median(),
                   decisive_correct_rate=usable.decisive_correct.mean(),
                   valid_four_model_comparisons=int(g.all_four_valid.sum()))
        for model in MODELS:
            row[f"winner_rate_{model}"] = g[f"winner_credit_{model}"].mean()
        srows.append(row)
    pd.DataFrame(srows).to_csv(out / "selection_summary.csv", index=False)
    diagnostic = pd.DataFrame(diagnostic_rows)
    diagnostic.to_csv(out / "third_order.csv", index=False)
    numeric = [c for c in diagnostic if c.startswith(("k3_", "r3_"))]
    if numeric:
        diagnostic.groupby("scenario_id")[numeric].agg(["mean", "std", "count"]).to_csv(
            out / "third_order_summary.csv")
    failure_counts = fits.loc[~fits.success, "error"].value_counts().to_dict()
    run_summary = {"datasets": len(records), "model_fits": len(fits),
                   "successful_fits": int(fits.success.sum()),
                   "valid_interval_sets": int(fits.interval_valid.sum()),
                   "failure_reasons": failure_counts,
                   "success_by_model": fits.groupby("model").success.agg(["sum", "count"]).to_dict("index"),
                   "max_truncation_loglik_change": fits.truncation_loglik_change.max()
                   if "truncation_loglik_change" in fits else None}
    write_json(out / "run_summary.json", run_summary)
    make_plots(out, summary)
    print(json.dumps(sanitize(run_summary), indent=2), flush=True)


def make_plots(out, summary):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    # Exportable scientific figures; no smoothing or invented observations.
    for parameter in ("mu", "r", "dep"):
        part = summary[summary.parameter == parameter]
        ids = list(part.scenario_id.unique())
        if not ids:
            continue
        height = max(4.0, len(ids) * .37 + 1)
        fig, ax = plt.subplots(figsize=(10.5, height), constrained_layout=True)
        colors = {"MPN0": "#718096", "MPN1": "#b84343", "MPNc": "#136f63", "CRTE": "#3766a2"}
        models = [m for m in MODELS if m in part.model.values]
        offsets = np.linspace(-.27, .27, len(models))
        for offset, model in zip(offsets, models):
            dat = part[part.model == model].set_index("scenario_id").reindex(ids)
            ax.plot(dat.rmse, np.arange(len(ids)) + offset, "o", ms=4,
                    color=colors[model], label=model)
        ax.set_yticks(np.arange(len(ids)), ids, fontsize=7)
        ax.set_xlabel(f"RMSE of {parameter} (successful fits only)")
        ax.set_ylabel("Scenario")
        ax.grid(axis="x", alpha=.2); ax.legend(frameon=False, ncol=len(models))
        ax.invert_yaxis()
        fig.savefig(out / f"rmse_{parameter}.png", dpi=180)
        fig.savefig(out / f"rmse_{parameter}.pdf")
        plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="command", required=True)
    run = sub.add_parser("run")
    run.add_argument("--suite", choices=("smoke", "core", "reverse", "small", "stress", "full"), default="smoke")
    run.add_argument("--reps", type=int)
    run.add_argument("--workers", type=int, default=min(8, os.cpu_count() or 1))
    run.add_argument("--seed", type=int, default=20260910)
    run.add_argument("--out", type=Path, required=True)
    run.add_argument("--starts", type=int, default=3)
    run.add_argument("--maxiter", type=int, default=250)
    run.add_argument("--no-intervals", action="store_true")
    report = sub.add_parser("report"); report.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    if args.command == "report":
        summarize(args.out)
    else:
        if not 1 <= args.workers <= 8:
            p.error("workers must lie between 1 and 8")
        reps = args.reps if args.reps is not None else (10 if args.suite == "smoke" else 500)
        if reps < 1 or args.starts < 1 or args.maxiter < 1:
            p.error("reps, starts and maxiter must be positive")
        cfg = FitConfig(starts=args.starts, maxiter=args.maxiter, intervals=not args.no_intervals)
        run_suite(args.suite, args.out, reps, args.workers, args.seed, cfg)


if __name__ == "__main__":
    main()
