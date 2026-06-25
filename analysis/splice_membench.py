#!/usr/bin/env python3
"""Splice freshly-measured bench cells into every figure-data JSON that holds
them (first-wins load_cells dedup means a stale copy in an earlier file would
otherwise win). Pass one or more source bench.json paths; matching (cell, scale)
entries replace the perf+RSS fields in-place across all figure data files."""
import json
import os
import sys

RES = "docs/manuscript/results"
FIGS = os.path.join(RES, "figures")

# Fields a fresh measurement supplies (wall / cpu / rss + their stats).
COPY = ["reps", "our_wall_ms", "up_wall_ms", "our_cpu_ms", "up_cpu_ms",
        "our_rss_mb", "up_rss_mb", "wall_ratio", "cpu_ratio", "rss_ratio",
        "our_samples", "up_samples", "our_wall_med", "our_wall_iqr",
        "up_wall_med", "up_wall_iqr", "our_cpu_med", "our_cpu_iqr",
        "up_cpu_med", "up_cpu_iqr", "our_rss_med", "our_rss_iqr",
        "up_rss_med", "up_rss_iqr", "wall_ratio_med", "cpu_ratio_med",
        "rss_ratio_med", "wall_ratio_ci_lo", "wall_ratio_ci_hi"]


def walk(x):
    if isinstance(x, dict):
        if "cell" in x:
            yield x
        for v in x.values():
            yield from walk(v)
    elif isinstance(x, list):
        for v in x:
            yield from walk(v)


# Build (cell, scale) -> fresh entry from the source bench files.
fresh = {}
for src in sys.argv[1:]:
    for e in walk(json.load(open(src))):
        if "cell" in e and "scale" in e:
            fresh[(e["cell"], e["scale"])] = e

# Every figure-data file load_cells reads.
files = [os.path.join(FIGS, "bench_multiscale.json"),
         os.path.join(RES, "full_validation", "bench.json"),
         os.path.join(FIGS, "bench_large_bamvcf.json"),
         os.path.join(FIGS, "bench_fastq.json"),
         os.path.join(FIGS, "bench_oom_large.json")]
lt = os.path.join(RES, "large_tier", "bench")
if os.path.isdir(lt):
    for g in sorted(os.listdir(lt)):
        files.append(os.path.join(lt, g, "bench.json"))

updated = 0
for fn in files:
    if not os.path.exists(fn):
        continue
    data = json.load(open(fn))
    changed = False
    for c in data if isinstance(data, list) else []:
        key = (c.get("cell"), c.get("scale"))
        if key in fresh:
            src = fresh[key]
            for f in COPY:
                if f in src:
                    c[f] = src[f]
            changed = True
            updated += 1
    if changed:
        json.dump(data, open(fn, "w"), indent=2)
        print("updated", os.path.relpath(fn, RES))
print("total entries updated:", updated)
