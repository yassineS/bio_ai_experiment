# Archived documentation

Historical snapshots of project state from the Phase 0/1 era (Oct 2025) and
the first big agent push (May 2026). Kept for reference; the live status
lives at the canonical locations linked below.

| Archived file | Original location | Replaced by |
|---|---|---|
| `PHASE1_COMPLETION.md` | repo root | [`tools/PORTING_STATUS.md`](../../tools/PORTING_STATUS.md) and [`README.md`](../../README.md) |
| `SESSION_SUMMARY_2025-10-21.md` | `tools/` | [`tools/PORTING_STATUS.md`](../../tools/PORTING_STATUS.md) |
| `TOOLS_IMPLEMENTATION_SUMMARY_2025-10-21.md` | `tools/IMPLEMENTATION_SUMMARY.md` | [`tools/PORTING_STATUS.md`](../../tools/PORTING_STATUS.md) |
| `GO_IMPLEMENTATION_SUMMARY.md` | `docs/` | [`tools/PORTING_STATUS.md`](../../tools/PORTING_STATUS.md) and [`CLAUDE.md`](../../CLAUDE.md) |

These files are **out of date**. They describe a state where 0-3 tools were
implemented; as of May 2026 the repo has 16+ tool ports with working
subsets. Don't trust the numbers in them.

A previous `SESSION_HANDOFF.md` doc lived in `docs/` but was removed in this
cleanup pass — the work it described has all landed on `main` and is now
captured by commit history and the canonical status files above.

A previous `top_50_packages_for_improvement.md` duplicate at the repo root
was also removed in this pass — the live version lives at
[`analysis/top_50_packages_for_improvement.md`](../../analysis/top_50_packages_for_improvement.md)
(itself superseded by
[`analysis/tool_ranking_2026.md`](../../analysis/tool_ranking_2026.md)).
