# bioformats (deprecated — see `pkg/htsgo`)

This directory is a transitional re-export layer for the htsgo
migration. The real format/index libraries live (or are landing) at
`pkg/htsgo/`; see [`pkg/htsgo/README.md`](../htsgo/README.md) for the
canonical inventory and [`docs/HTSGO_ROADMAP.md`](../../docs/HTSGO_ROADMAP.md)
for the migration plan.

## Where each format lives today

| Format          | Old import path (still works, deprecated) | New import path                              |
|-----------------|-------------------------------------------|-----------------------------------------------|
| `iohelper`      | `pkg/bioformats/iohelper`                 | `pkg/htsgo/iohelper`                         |
| `fasta`         | `pkg/bioformats/fasta`                    | `pkg/htsgo/fasta`                            |
| `fastq`         | `pkg/bioformats/fastq`                    | `pkg/htsgo/fastq`                            |
| `bed`           | `pkg/bioformats/bed`                      | `pkg/htsgo/bed`                              |
| `gff`           | `pkg/bioformats/gff`                      | `pkg/htsgo/gff`                              |
| `sam`           | `pkg/bioformats/sam`                      | `pkg/htsgo/sam`                              |
| `vcf`           | `pkg/bioformats/vcf`                      | `pkg/htsgo/vcf`                              |
| `bcf`           | `pkg/bioformats/bcf`                      | `pkg/htsgo/bcf`                              |

New code should import the `pkg/htsgo/` path directly. The shims here
will be deleted in the final PR of the migration (PR-I) by a single
mechanical rename-imports commit across the tree.

## Why a transitional layer

There are ~220 importers across `tools/`. Doing a single big-bang move
would inflate every migration PR with rename churn and make review
intractable. The shims let each move land on its own with `main`
green, deferring the rename sweep to one focused tail PR.

Each shim is a tiny `shim.go` file using Go type aliases
(`type Foo = htsgo.Foo`) and function variable re-exports
(`var NewReader = htsgo.NewReader`). No behaviour changes pass
through the shim — it's a pure forwarder.
