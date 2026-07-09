# Porting status — where to look

This file used to hold a hand-maintained per-subcommand feature inventory.
That inventory drifted out of date faster than it could be kept honest, so the
per-tool / per-subcommand status now lives in two source-of-truth documents,
and this file is a thin redirect to them. (Several docs — `CLAUDE.md`, the
docs map, `PROJECT_STATUS.md`, `pkg/htsgo/README.md`, `docs/HTSGO_ROADMAP.md`,
and per-tool progress notes — link here, so it is kept as a stable pointer
rather than deleted.)

## Source of truth

| You want… | Look in |
| --- | --- |
| Top-level per-tool completion table (the summary view) | [`../PROJECT_STATUS.md`](../PROJECT_STATUS.md) |
| Authoritative per-tool / per-subcommand **parity gap** detail | [`../docs/PARITY_ROADMAP.md`](../docs/PARITY_ROADMAP.md) |
| Per-tool usage + parity notes | `tools/<tool>/README.md` |
| htsgo format-library parity matrix | [`../pkg/htsgo/README.md`](../pkg/htsgo/README.md) and [`../docs/HTSGO_ROADMAP.md`](../docs/HTSGO_ROADMAP.md) |
| The docs map (which document owns which kind of status) | [`../docs/README.md`](../docs/README.md) |

## Per-subcommand inventories

The rich per-subcommand feature lists (e.g. every `samtools` / `bcftools`
subcommand and which flags are ported) are maintained **inline in
`docs/PARITY_ROADMAP.md`**, under each tool's section, alongside the remaining
gaps. That co-location keeps the inventory and the gap list from disagreeing —
the failure mode that retired the standalone table here.

## Note on repository structure

This is a **single Go module** (`github.com/yassineS/bio_ai_experiment`); there
is no per-tool `go.mod`. Tool logic lives in `tools/<tool>/pkg/<tool>/` with the
CLI entry point in `tools/<tool>/cmd/<tool>/main.go`, and tests sit next to the
code as `*_test.go`. See `CLAUDE.md` for the full conventions.
