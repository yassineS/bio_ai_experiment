# Project status

**Last refreshed:** 2026-05-14

This file used to describe the Phase 0/1 setup-and-analysis period (Oct 2025)
and said "0 tools implemented", which has been wrong for some time. The
historical text has been moved to
[`docs/archive/PHASE1_COMPLETION.md`](docs/archive/PHASE1_COMPLETION.md). This
file is now a thin pointer to the live status.

## Current snapshot

- **16 tool ports** with a working subset of features. Per-tool coverage and
  feature lists live in [`tools/PORTING_STATUS.md`](tools/PORTING_STATUS.md).
- Single Go module at the repo root (`github.com/yassineS/bio_ai_experiment`),
  **zero third-party Go dependencies**. See [`CLAUDE.md`](CLAUDE.md) for the
  repo layout and conventions.
- POSIX-compliant CLIs are a parity requirement; see
  [`docs/CLI_CONVENTIONS.md`](docs/CLI_CONVENTIONS.md).
- CI workflow is currently **disabled** (manual-only via `workflow_dispatch`)
  while the project iterates heavily; agents run `gofmt`/`go vet`/`go build`/
  `go test -race -cover` + `markdownlint` locally and document the output in
  each PR description.

## Where to look next

| Question | Doc |
|---|---|
| What's the per-tool state? | [`tools/PORTING_STATUS.md`](tools/PORTING_STATUS.md) |
| Which tools should we port next? | [`analysis/tool_ranking_2026.md`](analysis/tool_ranking_2026.md) |
| How is the repo organised? | [`CLAUDE.md`](CLAUDE.md) |
| What are the CLI rules? | [`docs/CLI_CONVENTIONS.md`](docs/CLI_CONVENTIONS.md) |
| How do I use tool X? | `tools/<tool>/README.md` |
| What's the project's overall pitch? | [`README.md`](README.md) |
