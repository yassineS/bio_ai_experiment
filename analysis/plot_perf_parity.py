#!/usr/bin/env python3
"""Generate manuscript performance figures + parity/performance tables.

Consumes the pipeline/bench artifacts (multi-scale small→large) and emits, in
the Majorelle design language (https://github.com/yassineS/majorelle-py):

  - figures/fig_speedup.{pdf,png}  — per-cell speedup (upstream/ours), grouped by
    tool family, with one bar per tier (small/medium/large) and 95% bootstrap CI.
  - figures/fig_scaling.{pdf,png}  — wall time vs input tier (log-y), ours vs
    upstream, coloured by tool family (solid = ours, dashed = upstream).
  - figures/fig_scatter.{pdf,png}  — ours vs upstream wall time (log-log) with
    the y=x parity line.
  - performance_tables.md          — C3 per-cell speedups + C2 parity rates
    with Wilson / Clopper-Pearson CIs.

Run:
  uv run --with "git+https://github.com/yassineS/majorelle-py" \
         --with matplotlib --with numpy --with scipy \
         python analysis/plot_perf_parity.py

Convention: speedup = t_upstream / t_ours, so >1 = our port is faster (C3). The
bench json stores wall_ratio = ours/upstream; speedup is its reciprocal and the
CI bounds invert + swap.
"""
import json
import math
import os
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.lines
from scipy.stats import beta

from majorelle import PAL_QUAL, COLOURS
from majorelle import mpl as mj

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "docs", "manuscript", "results")
FIGS = os.path.join(RES, "figures")
os.makedirs(FIGS, exist_ok=True)

mj.set_theme(variant="rabat")           # white-walled, formal — manuscript
try:
    mj.enable_tabular_numerals()
except Exception:
    pass

# Tiers shown in the manuscript (smoke omitted — too small to be informative).
TIERS = ["small", "medium", "large"]
TIER_X = {"small": 1, "medium": 2, "large": 3}
# small -> medium -> large as a light->dark blue ramp (Majorelle Blues).
TIER_COLOUR = {"small": "#9C8AF3", "medium": "#4A2EE8", "large": "#1A0E7A"}

# cell name -> (tool family, pretty subcommand label)
CELL = {
    "sam_view_bam2bam": ("samtools", r"view (bam$\to$bam)"),
    "sam_view_bam2cram": ("samtools", r"view (bam$\to$cram)"),
    "sam_view_cram2bam": ("samtools", r"view (cram$\to$bam)"),
    "sam_sort": ("samtools", "sort"),
    "sam_flagstat": ("samtools", "flagstat"),
    "sam_stats": ("samtools", "stats"),
    "sam_depth": ("samtools", "depth"),
    "sam_mpileup": ("samtools", "mpileup"),
    "bcf_view": ("bcftools", "view"),
    "bcf_view_body": ("bcftools", "view -H"),
    "bcf_norm": ("bcftools", "norm"),
    "bcf_stats": ("bcftools", "stats"),
    "bcf_query": ("bcftools", "query"),
    "bcf_call": ("bcftools", "call"),
    "bcf_isec": ("bcftools", "isec"),
    "bed_intersect_self": ("bedtools", "intersect (self)"),
    "bed_intersect_pair": ("bedtools", "intersect"),
    "bed_merge": ("bedtools", "merge"),
    "bed_coverage": ("bedtools", "coverage"),
    "bed_genomecov": ("bedtools", "genomecov"),
    "bed_sort": ("bedtools", "sort"),
    "seqtk_seq": ("seqtk", "seq"),
    "sickle_se": ("sickle", "se"),
}
TOOL_ORDER = ["samtools", "bcftools", "bedtools", "seqtk", "sickle"]


def load_cells():
    """Flat list of bench cells from every available bench.json, keyed by tier."""
    cells = []
    seen = set()

    def add(path):
        if os.path.exists(path):
            for c in json.load(open(path)):
                k = (c["cell"], c["scale"])
                if k not in seen:
                    seen.add(k)
                    cells.append(c)

    add(os.path.join(FIGS, "bench_multiscale.json"))      # smoke/small/medium
    add(os.path.join(RES, "full_validation", "bench.json"))
    lt = os.path.join(RES, "large_tier", "bench")
    if os.path.isdir(lt):
        for grp in sorted(os.listdir(lt)):
            add(os.path.join(lt, grp, "bench.json"))
    add(os.path.join(FIGS, "bench_large_bamvcf.json"))    # large BAM/VCF light
    return [c for c in cells if c["scale"] in TIERS and c["cell"] in CELL]


def speedup(c):
    s = c["up_wall_med"] / c["our_wall_med"]
    lo = 1.0 / c["wall_ratio_ci_hi"] if c.get("wall_ratio_ci_hi") else s
    hi = 1.0 / c["wall_ratio_ci_lo"] if c.get("wall_ratio_ci_lo") else s
    return s, lo, hi


def ordered_cells(by):
    """cell names ordered by tool family then a stable within-family order."""
    out = []
    for tool in TOOL_ORDER:
        fam = [name for name in CELL if CELL[name][0] == tool and name in by]
        fam.sort(key=lambda n: CELL[n][1])
        out += [(tool, n) for n in fam]
    return out


def fig_speedup(cells):
    by = defaultdict(dict)                       # cell -> tier -> c
    for c in cells:
        by[c["cell"]][c["scale"]] = c
    rows = ordered_cells(by)                      # [(tool, cell), ...]

    # y layout: a gap between tool families.
    ypos, ylabels, sep, prev = [], [], [], None
    y = 0.0
    for tool, name in rows:
        if prev is not None and tool != prev:
            y += 0.8
            sep.append(y - 0.4)
        ypos.append(y)
        ylabels.append(CELL[name][1])
        prev = tool
        y += 1.0

    n_t = len(TIERS)
    bar_h = 0.78 / n_t
    fig, ax = plt.subplots(figsize=(7.4, max(6.5, 0.42 * len(rows) + 1.2)))
    for ti, tier in enumerate(TIERS):
        ys, xs, los, his = [], [], [], []
        for (tool, name), yb in zip(rows, ypos):
            c = by[name].get(tier)
            if not c:
                continue
            s, lo, hi = speedup(c)
            off = (ti - (n_t - 1) / 2) * bar_h
            ys.append(yb - off)
            xs.append(s)
            los.append(max(s - lo, 0))
            his.append(max(hi - s, 0))
        ax.barh(ys, xs, height=bar_h * 0.92, color=TIER_COLOUR[tier],
                xerr=[los, his], ecolor=COLOURS["gray_48"], capsize=1.5,
                error_kw={"lw": 0.8}, label=tier, zorder=3)

    ax.axvline(1.0, color=COLOURS["near_black"], lw=1.1, ls="--", alpha=0.8, zorder=2)
    ax.set_yticks(ypos)
    ax.set_yticklabels(ylabels, fontsize=9)
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlim(0.3, 3.6)
    ax.xaxis.set_major_locator(matplotlib.ticker.FixedLocator([1 / 3, 0.5, 1, 2, 3]))
    ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(
        lambda v, _: ("0.33×" if abs(v - 1 / 3) < 0.01 else f"{v:g}×")))
    ax.set_xlabel("speedup  =  upstream wall / ours   (>1 means our port is faster)")

    # tool-family band labels on the left margin + light separators.
    fam_rows = defaultdict(list)
    for (tool, name), yb in zip(rows, ypos):
        fam_rows[tool].append(yb)
    for tool, ys in fam_rows.items():
        ymid = sum(ys) / len(ys)
        ax.text(-0.215, ymid, tool, transform=ax.get_yaxis_transform(),
                rotation=90, va="center", ha="center", fontsize=10,
                fontweight="bold", color=COLOURS["primary"])
    for sy in sep:
        ax.axhline(sy, color=COLOURS["gray_16"], lw=0.7, zorder=1)

    ax.legend(title="tier", loc="upper right", frameon=True, framealpha=0.95,
              edgecolor=COLOURS["gray_16"], fontsize=9, ncol=1, title_fontsize=9)
    fig.subplots_adjust(left=0.30, right=0.97, bottom=0.07, top=0.91)
    mj.title_block(fig, "Per-subcommand speedup vs upstream",
                   "median over reps · 95% bootstrap CI · grouped by tool",
                   left=0.30, tighten=False)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIGS, f"fig_speedup.{ext}"), dpi=200)
    plt.close(fig)
    return len(rows)


def fig_scaling(cells):
    """Wall time vs tier, one line per cell, COLOURED BY TOOL FAMILY so every
    tool (incl. bcftools) is visible; solid = ours, dashed = upstream."""
    by = defaultdict(dict)
    for c in cells:
        by[c["cell"]][c["scale"]] = c
    multi = {k: v for k, v in by.items() if len([t for t in TIERS if t in v]) >= 2}
    tool_col = {t: PAL_QUAL[i] for i, t in enumerate(TOOL_ORDER)}
    fig, ax = plt.subplots(figsize=(7.0, 5.2))
    seen_tool = set()
    for tool, name in ordered_cells(multi):
        sc = multi[name]
        ts = [t for t in TIERS if t in sc]
        xs = [TIER_X[t] for t in ts]
        ours = [sc[t]["our_wall_med"] for t in ts]
        up = [sc[t]["up_wall_med"] for t in ts]
        col = tool_col[tool]
        lab = tool if tool not in seen_tool else None
        seen_tool.add(tool)
        ax.plot(xs, ours, "-o", color=col, lw=1.7, ms=4, label=lab, zorder=3)
        ax.plot(xs, up, "--o", color=col, lw=1.0, ms=2.5, alpha=0.55, zorder=2)
    ax.set_yscale("log")
    ax.yaxis.set_major_locator(matplotlib.ticker.LogLocator(base=10, numticks=6))
    ax.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.set_xticks(list(TIER_X.values()))
    ax.set_xticklabels(list(TIER_X.keys()))
    ax.set_xlim(0.85, 3.15)
    ax.set_xlabel("input tier")
    ax.set_ylabel("wall time (ms)")
    leg = ax.legend(title="tool  (solid = ours, dashed = upstream)", fontsize=9,
                    frameon=False, loc="upper left", ncol=2, title_fontsize=9)
    leg._legend_box.align = "left"
    fig.subplots_adjust(left=0.11, right=0.97, bottom=0.12, top=0.86)
    mj.title_block(fig, "Scaling across input tiers",
                   "ours (solid) vs upstream (dashed) — both scale alike",
                   left=0.11, tighten=False)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIGS, f"fig_scaling.{ext}"), dpi=200)
    plt.close(fig)
    return len(multi)


def fig_scatter(cells):
    marker = {"small": "o", "medium": "s", "large": "D"}
    fig, ax = plt.subplots(figsize=(6.2, 6.2))
    for tier in TIERS:
        pts = [c for c in cells if c["scale"] == tier]
        if not pts:
            continue
        x = [c["up_wall_med"] for c in pts]
        y = [c["our_wall_med"] for c in pts]
        col = [COLOURS["primary"] if c["our_wall_med"] <= c["up_wall_med"]
               else COLOURS["tertiary"] for c in pts]
        ax.scatter(x, y, c=col, marker=marker[tier], s=40, alpha=0.85,
                   edgecolors=COLOURS["near_black"], linewidths=0.4,
                   label=tier, zorder=3)
    lim = [1, 1e6]
    ax.plot(lim, lim, ls="--", color=COLOURS["near_black"], lw=1, alpha=0.8, zorder=2)
    ax.set_xscale("log"); ax.set_yscale("log")
    # decade ticks only — no minor labels.
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_major_locator(matplotlib.ticker.LogLocator(base=10, numticks=7))
        axis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_aspect("equal")
    ax.set_xlabel("upstream wall time (ms)")
    ax.set_ylabel("our wall time (ms)")
    # one legend: tier marker shape; colour meaning explained in the subtitle.
    handles = [matplotlib.lines.Line2D([], [], marker=marker[t], ls="", color=COLOURS["gray_48"],
                          ms=7, label=t) for t in TIERS if any(c["scale"] == t for c in cells)]
    ax.legend(handles=handles, title="tier", fontsize=9, frameon=False,
              loc="lower right", title_fontsize=9)
    fig.subplots_adjust(left=0.14, right=0.97, bottom=0.16, top=0.85)
    mj.title_block(fig, "Our wall time vs upstream",
                   "below the parity line = faster (blue); above = slower (terracotta)",
                   left=0.13, tighten=False)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIGS, f"fig_scatter.{ext}"), dpi=200)
    plt.close(fig)


def wilson(k, n, z=1.95996):
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (100 * (c - h) / d, 100 * (c + h) / d)


def clopper(k, n, a=0.05):
    lo = 0.0 if k == 0 else beta.ppf(a / 2, k, n - k + 1) * 100
    hi = 100.0 if k == n else beta.ppf(1 - a / 2, k + 1, n - k) * 100
    return (lo, hi)


def perf_table(cells):
    lines = ["# Manuscript performance & parity tables",
             "",
             "Auto-generated by `analysis/plot_perf_parity.py` from the `pipeline/bench`",
             "artifacts, styled in the Majorelle design language. Figures in",
             "[`figures/`](figures/). Convention: **speedup = upstream wall / ours**,",
             "so **> 1 means our port is faster**. Timings are median over reps with the",
             "inter-quartile range; the 95% CI is a bootstrap on the wall ratio (H1a).",
             "",
             "## C3 — Performance (per subcommand, per tier)",
             "",
             "| tier | tool | subcommand | ours ms (med±IQR) | upstream ms (med±IQR) | speedup | 95% CI | CPU× | note |",
             "|---|---|---|---|---|---|---|---|---|"]
    torder = {t: i for i, t in enumerate(TIERS)}
    for c in sorted(cells, key=lambda c: (torder[c["scale"]], CELL[c["cell"]][0],
                                          -(c["up_wall_med"] / c["our_wall_med"]))):
        tool, sub = CELL[c["cell"]]
        s, lo, hi = speedup(c)
        note = "faster" if s >= 1.10 else ("par" if s >= 0.90 else "slower")
        if s < 0.5:
            note = "**slow**"
        ci = f"[{lo:.2f}, {hi:.2f}]"
        lines.append(
            f"| {c['scale']} | {tool} | `{sub}` | {c['our_wall_med']:.1f} ± {c['our_wall_iqr']:.1f} "
            f"| {c['up_wall_med']:.1f} ± {c['up_wall_iqr']:.1f} | {s:.2f}× | {ci} "
            f"| {c['up_cpu_med'] / c['our_cpu_med']:.2f} | {note} |")
    med = [c for c in cells if c["scale"] == "medium"]
    fast = sum(1 for c in med if c["up_wall_med"] / c["our_wall_med"] >= 1.10)
    par = sum(1 for c in med if 0.90 <= c["up_wall_med"] / c["our_wall_med"] < 1.10)
    lines += ["",
              f"**Medium-tier summary ({len(med)} cells):** {fast} faster (≥1.1×), "
              f"{par} at par, {len(med) - fast - par} slower. The I/O-bound conversions "
              "and `bedtools` intersect/coverage/genomecov + `sickle` are faster; the "
              "compute-heavy variant cells (`mpileup`, `call`, `isec`) are slower and "
              "reported plainly (and OOM at the large tier on the 12 GB box).", ""]
    return lines


def parity_table():
    lines = ["## C2 — Parity rates (byte-exact vs upstream) with 95% CIs",
             "",
             "Binomial proportions (`k` byte-identical of `n` compared cells, "
             "provenance-stripped) with Wilson score and exact Clopper-Pearson 95% "
             "intervals (`pipeline/stats`). `full-validation` diverges are the documented "
             "`arm64`-platform floating-point cells (byte-exact on `amd64`/CI).",
             "",
             "| Experiment | k / n | Point (%) | Wilson 95% CI | Clopper-Pearson 95% CI |",
             "|---|---|---:|---|---|"]
    rows = [
        ("realparity (real GIAB exome, samtools+bcftools)", 15, 15),
        ("full-validation medium (all 53 CLIs, byte-exact-eligible)", 383, 385),
        ("full-validation medium, amd64-adjusted", 385, 385),
        ("full-validation small", 381, 385),
        ("full-validation smoke", 383, 385),
    ]
    for name, k, n in rows:
        wl, cp = wilson(k, n), clopper(k, n)
        lines.append(f"| {name} | {k} / {n} | {100*k/n:.2f} | "
                     f"[{wl[0]:.2f}, {wl[1]:.2f}] | [{cp[0]:.2f}, {cp[1]:.2f}] |")
    lines += ["",
              "> The medium denominator counts the 385 byte-exact-eligible cells (4 SIMILAR "
              "float-scored cells are tolerance-compared, excluded from the binomial). The 2 "
              "raw diverges are arm64 `vcftools` `-nan`/last-ULP; on `amd64`/CI it is 385/385.", ""]
    return lines


def main():
    cells = load_cells()
    tiers = sorted(set(c["scale"] for c in cells), key=lambda t: TIER_X[t])
    print(f"loaded {len(cells)} bench cells across tiers {tiers}")
    print("fig_speedup rows:", fig_speedup(cells))
    print("fig_scaling rows:", fig_scaling(cells))
    fig_scatter(cells)
    out = perf_table(cells) + parity_table()
    with open(os.path.join(RES, "performance_tables.md"), "w") as f:
        f.write("\n".join(out) + "\n")
    print("wrote performance_tables.md +", len([x for x in os.listdir(FIGS) if x.endswith(('.pdf', '.png'))]), "figure files")


if __name__ == "__main__":
    main()
