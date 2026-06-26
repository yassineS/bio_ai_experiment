# Documentation map

This file is the index of the project's documentation. Its main job is to
state **which document owns which kind of information**, so status and
progress content lives in exactly one place and contributors don't re-fork it.

## Single source of truth — who owns what

| You want… | Go to | Notes |
|-----------|-------|-------|
| Top-level "distance to done" completion table | [`../PROJECT_STATUS.md`](../PROJECT_STATUS.md) | Skimmable per-tool % + biggest boulders. **Owns** the summary view. |
| Authoritative per-tool parity gap list | [`PARITY_ROADMAP.md`](PARITY_ROADMAP.md) | The canonical, detailed "what's missing per subcommand/flag" tracker. **Owns** the gap detail. |
| Per-subcommand feature inventory + test/coverage notes | [`PARITY_ROADMAP.md`](PARITY_ROADMAP.md) + per-tool [`../tools/<tool>/README.md`](../tools/) | Feature-by-feature gaps live in the parity roadmap; per-tool READMEs carry the usage + parity notes. |
| Validated-parity test methodology + skip ledger | [`FULL_VALIDATION.md`](FULL_VALIDATION.md) + [`CONFORMANCE.md`](CONFORMANCE.md) | How parity is proven against upstream corpora, and the validation run recipes. |
| Upstream bugs we fixed on port | [`UPSTREAM_BUGS.md`](UPSTREAM_BUGS.md) | Deviations from upstream that are intentional fixes. |
| Documented intentional CLI differences | [`CLI_CONVENTIONS.md`](CLI_CONVENTIONS.md) + per-tool [`../tools/<tool>/README.md`](../tools/) | The canonical CLI flag spec; per-tool READMEs note where flags/behaviour intentionally differ. |
| Which **new** tools to port next | [`../analysis/tool_ranking_2026.md`](../analysis/tool_ranking_2026.md) | A priority list for *new* ports, not a deprioritise filter for existing ones. |
| How to use a specific tool | `../tools/<tool>/README.md` | Per-tool usage + parity notes. |
| Repo layout, build/test commands, deps policy | [`../CLAUDE.md`](../CLAUDE.md) | The orientation doc for contributors and AI agents. |
| How to contribute (PRs, issues, setup) | [`../CONTRIBUTING.md`](../CONTRIBUTING.md) | Workflow and coding standards. |

If two documents disagree on status, **`PROJECT_STATUS.md` (summary) and
`PARITY_ROADMAP.md` (detail) win.** Everything else links to them rather than
restating numbers.

## Design and reference docs (own their topic, not status)

- [`GOLANG_GUIDE.md`](GOLANG_GUIDE.md) — Go patterns and best practices used here.
- [`CLI_CONVENTIONS.md`](CLI_CONVENTIONS.md) — the canonical CLI flag spec
  (POSIX short + GNU long via `pkg/cliflag`).
- [`CRAM_DESIGN.md`](CRAM_DESIGN.md) / [`CRAM_ROADMAP.md`](CRAM_ROADMAP.md) —
  the CRAM port's up-front decisions and remaining work.
- [`HTSGO_ROADMAP.md`](HTSGO_ROADMAP.md) — consolidating format/index code
  into the shared `pkg/htsgo` htslib-equivalent library.
- [`PLUGIN_PROTOCOL.md`](PLUGIN_PROTOCOL.md) — the bcftools plugin subprocess
  protocol.
- [`../pkg/htsgo/README.md`](../pkg/htsgo/README.md) — format library docs.
- [`METRICS.md`](METRICS.md) — scope, lines-of-code, and speed vs the originals
  (the manuscript metrics); the benchmark harness lives in
  [`../pipeline/bench`](../pipeline/bench).

## Agent role descriptions

The AI agents that build this repo are described under
[`../.github/agents/`](../.github/agents/): tool-analysis, golang-recoding,
testing, and documentation roles, plus a
[coordination guide](../.github/agents/COORDINATION.md). The mcp-server role is
**deprecated** — MCP servers are not being built (the project ships drop-in
POSIX CLIs instead); see
[`../.github/agents/mcp-server-agent.md`](../.github/agents/mcp-server-agent.md).

## Documentation conventions

- Keep Markdown well-formed: the markdown lint over `**/*.md` is part of the
  (currently `workflow_dispatch`-only) CI job and is run locally on PRs.
- Put status/progress facts in the owning document above and **link** to them
  elsewhere instead of copying — duplicated status is how docs go stale.
- Document exported Go identifiers with complete-sentence doc comments
  (see [`GOLANG_GUIDE.md`](GOLANG_GUIDE.md)).
