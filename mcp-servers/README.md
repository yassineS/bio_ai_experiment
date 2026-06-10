# MCP Servers — DESCOPED

**This initiative has been descoped. No MCP (Model Context Protocol) servers
are being built for this project, and none are planned.**

The project originally intended to ship a per-tool MCP server so the recoded
bioinformatics tools could be driven by Large Language Models. The owner has
since decided **not** to pursue that direction. This directory is retained only
as a marker; do not start building MCP servers here.

## What the interface is instead

Each recoded tool is exposed as a **standalone, POSIX-compliant command-line
program** that is **retro-compatible (a drop-in replacement) with its upstream
tool**. That means:

- POSIX short options (`-i`, `-q`, …) and GNU-style long options (`--input`,
  `--quality`, …), with getopt-style flag bundling.
- `--` to end option parsing, `-` for stdin/stdout, and clean exit codes.
- The upstream tool's flag surface preserved wherever practical, so existing
  pipelines and scripts work unchanged.

This CLI work is largely complete — `cliflag.Parse` (getopt bundling) plus
per-tool flag retro-compatibility has rolled out across all tools.

## Where to look

- [`../docs/CLI_CONVENTIONS.md`](../docs/CLI_CONVENTIONS.md) — the canonical CLI
  flag and behavior spec.
- [`../PROJECT_STATUS.md`](../PROJECT_STATUS.md) — per-tool completion status.
