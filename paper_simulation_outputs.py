#!/usr/bin/env python3
"""Regenerate the manuscript's simulation tables and figures from a run folder.

Reads the outputs written by ``softclosure.py run``/``report`` (``fits.csv``,
``parameter_summary.csv`` and ``model_selection.csv``) and writes, into
``OUT/paper/``:

* ``fit_status.csv``                 - Table "Fit success and interval availability"
* ``paired_selection_pooled.csv``    - pooled paired AIC credit by block and T
* ``four_candidate_selection.csv``   - Table "Pooled four-candidate AIC selection credit"
* ``misspecification_cell.csv``      - Table "Correct and misspecified partial-closure fits"
* ``full_misspecification.{png,pdf}``  - Figure: mean estimates in that cell
* ``full_model_selection.{png,pdf}``   - Figure: paired AIC discrimination

Pooled credits weight scenarios by their numbers of valid comparisons, i.e. they
average replicate-level credits, exactly as described in the manuscript.
Nothing here re-fits a model; every number comes from the run folder.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

MODELS = ("MPN0", "MPN1", "MPNc", "CRTE")
PAPER_CELL = dict(S=200, T=7, true_mu=2.0, true_r=0.5, true_dep=0.8)


def fit_status(fits):
    rows = []
    for model in MODELS:
        g = fits[fits.model == model]
        rows.append(dict(model=model, attempted=len(g), successful=int(g.success.sum()),
                         success_percent=100 * g.success.mean() if len(g) else np.nan,
                         interval_sets=int(g.interval_valid.sum())))
    total = dict(model="All", attempted=len(fits), successful=int(fits.success.sum()),
                 success_percent=100 * fits.success.mean(),
                 interval_sets=int(fits.interval_valid.sum()))
    return pd.DataFrame(rows + [total])


def paired_selection(selection):
    usable = selection[selection.pair_valid & selection.mechanisms_distinct]
    rows = []
    for (block, T), g in usable.groupby(["block", "T"], sort=True):
        per = g.groupby("scenario_id").correct_selection_credit.mean()
        credit = g.correct_selection_credit.to_numpy(float)
        rows.append(dict(block=block, T=T, scenarios=len(per), valid_pairs=len(g),
                         pooled_credit_percent=100 * credit.mean(),
                         pooled_mcse_percent=100 * credit.std(ddof=1) / np.sqrt(len(credit))
                         if len(credit) > 1 else np.nan,
                         scenario_min_percent=100 * per.min(),
                         scenario_max_percent=100 * per.max()))
    return pd.DataFrame(rows)


def four_candidate(selection):
    valid = selection[selection.all_four_valid.astype(bool)]
    rows = []
    for (block, T), g in valid.groupby(["block", "T"], sort=True):
        row = dict(block=block, T=T, valid=len(g))
        for model in MODELS:
            row[f"{model}_percent"] = 100 * g[f"winner_credit_{model}"].mean()
        rows.append(row)
    return pd.DataFrame(rows)


def cell(summary, cell_spec, generator):
    mask = summary.dgm.eq(generator)
    for key, value in cell_spec.items():
        mask &= np.isclose(summary[key].astype(float), value)
    return summary[mask & summary.model.isin(("MPNc", "CRTE"))]


def misspecification_table(summary, cell_spec):
    rows = []
    for generator in ("MPNc", "CRTE"):
        part = cell(summary, cell_spec, generator)
        for model in ("MPNc", "CRTE"):
            sub = part[part.model == model].set_index("parameter")
            if sub.empty:
                continue
            rows.append(dict(generator=generator, fitted_model=model,
                             mean_abundance=sub.loc["mu", "mean"],
                             mean_r=sub.loc["r", "mean"],
                             mean_dependence=sub.loc["dep", "mean"],
                             successful=int(sub.loc["mu", "successful_fits"]),
                             abundance_relative_bias_percent=sub.loc["mu", "relative_bias_percent"],
                             r_relative_bias_percent=sub.loc["r", "relative_bias_percent"]))
    return pd.DataFrame(rows)


def plot_misspecification(summary, cell_spec, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    combos = [(g, m) for g in ("MPNc", "CRTE") for m in ("MPNc", "CRTE")]
    params = [("mu", "Visit-level abundance"), ("r", "Per-pass detection r"),
              ("dep", "Dependence (c or $\\phi$)")]
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.6), constrained_layout=True)
    colors = {"MPNc": "#136f63", "CRTE": "#3766a2"}
    found = False
    for ax, (param, label) in zip(axes, params):
        for i, (generator, model) in enumerate(combos):
            row = cell(summary, cell_spec, generator)
            row = row[(row.model == model) & (row.parameter == param)]
            if row.empty:
                continue
            found = True
            mean, mcse = float(row["mean"].iloc[0]), float(row["bias_mcse"].iloc[0])
            ax.errorbar(i, mean, yerr=1.96 * mcse, fmt="o", color=colors[model],
                        capsize=3, ms=5)
        target = {"mu": cell_spec["true_mu"], "r": cell_spec["true_r"],
                  "dep": cell_spec["true_dep"]}[param]
        ax.axhline(target, ls="--", lw=1, color="0.4")
        ax.set_xticks(range(4), [f"gen {g}\nfit {m}" for g, m in combos], fontsize=8)
        ax.set_title(label, fontsize=10)
        ax.grid(axis="y", alpha=.2)
    if not found:
        plt.close(fig)
        return False
    fig.savefig(out / "full_misspecification.png", dpi=200)
    fig.savefig(out / "full_misspecification.pdf")
    plt.close(fig)
    return True


def plot_selection(selection, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    usable = selection[selection.pair_valid & selection.mechanisms_distinct
                       & selection.block.isin(("core", "reverse"))]
    if usable.empty:
        return False
    agg = (usable.groupby(["dgm", "T", "true_mu", "true_r", "true_dep"])
           .correct_selection_credit.agg(["mean", "std", "count"]).reset_index())
    agg["mcse"] = agg["std"] / np.sqrt(agg["count"])
    gens, visits = ("MPNc", "CRTE"), sorted(agg["T"].unique())
    fig, axes = plt.subplots(len(gens), len(visits), figsize=(4.2 * len(visits), 6.4),
                             sharey=True, squeeze=False, constrained_layout=True)
    styles = {}
    palette = ["#136f63", "#3766a2", "#b84343", "#a0761f"]
    for (mu, r) in sorted(set(zip(agg.true_mu, agg.true_r))):
        styles[(mu, r)] = palette[len(styles) % len(palette)]
    for i, generator in enumerate(gens):
        for k, T in enumerate(visits):
            ax = axes[i, k]
            part = agg[(agg.dgm == generator) & (agg["T"] == T)]
            for (mu, r), g in part.groupby(["true_mu", "true_r"]):
                g = g.sort_values("true_dep")
                lo = np.clip(g["mean"] - 1.96 * g.mcse, 0, 1)
                hi = np.clip(g["mean"] + 1.96 * g.mcse, 0, 1)
                ax.errorbar(g.true_dep, g["mean"], yerr=[g["mean"] - lo, hi - g["mean"]],
                            marker="o", ms=4, capsize=2, color=styles[(mu, r)],
                            label=f"$\\mu$={mu:g}, r={r:g}")
            ax.set_title(f"{generator} generator, T={T}", fontsize=10)
            ax.set_ylim(0.4, 1.02)
            ax.axhline(0.5, ls=":", lw=1, color="0.5")
            ax.grid(alpha=.2)
            if i == len(gens) - 1:
                ax.set_xlabel("Generating dependence (c or $\\phi$)")
            if k == 0:
                ax.set_ylabel("Correct-selection credit")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="outside lower center", ncol=len(handles), frameon=False)
    fig.savefig(out / "full_model_selection.png", dpi=200)
    fig.savefig(out / "full_model_selection.pdf")
    plt.close(fig)
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True, help="softclosure.py run folder")
    args = parser.parse_args()
    run = args.out
    fits = pd.read_csv(run / "fits.csv")
    summary = pd.read_csv(run / "parameter_summary.csv")
    selection = pd.read_csv(run / "model_selection.csv")
    paper = run / "paper"
    paper.mkdir(exist_ok=True)
    tables = {"fit_status": fit_status(fits),
              "paired_selection_pooled": paired_selection(selection),
              "four_candidate_selection": four_candidate(selection),
              "misspecification_cell": misspecification_table(summary, PAPER_CELL)}
    for name, table in tables.items():
        table.to_csv(paper / f"{name}.csv", index=False)
        print(f"\n== {name} ==\n{table.round(3).to_string(index=False) if len(table) else '(no rows in this run)'}")
    print("\nfigure full_misspecification:",
          "written" if plot_misspecification(summary, PAPER_CELL, paper) else "cell not in this run")
    print("figure full_model_selection:",
          "written" if plot_selection(selection, paper) else "no core/reverse scenarios in this run")


if __name__ == "__main__":
    main()
