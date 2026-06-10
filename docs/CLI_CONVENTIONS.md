# CLI Design Conventions for Bio AI Experiment Tools

This document outlines the command-line interface (CLI) design conventions used across all tools in the Bio AI Experiment repository. These conventions ensure consistency and usability across all reimplemented bioinformatics tools.

## Project goal: POSIX-compliant CLIs at full parity

A tool in this repository is **only considered "complete"** once it has both:

1. **100% feature parity** with the original tool it ports; and
2. A **POSIX-compliant CLI** as specified in this document — POSIX-style short
   options, GNU-style long options, `-` for stdin/stdout, `--` to end option
   parsing, predictable exit codes, and `-h/--help` + `-v/--version`.

Until both conditions are met, the tool is listed as **Partial** in
[`tools/PORTING_STATUS.md`](../tools/PORTING_STATUS.md). Where a POSIX-compliant
flag conflicts with an upstream short-form alias, the POSIX form wins; the
upstream alias may still be accepted as a long option (e.g. `--qual-type`).

## General Principles

1. **POSIX-style flags**: Use both short (single dash, single letter) and long (double dash, full word) options
2. **Consistency**: Use the same option names for similar functionality across all tools
3. **Discoverability**: Provide helpful error messages and comprehensive help text
4. **Backward compatibility**: Where possible, maintain compatibility with original tool options, but never at the cost of POSIX compliance once a tool is at parity
5. **`--` ends option parsing** and `-` (alone) means stdin/stdout; both are required for compliance
6. **POSIX short-flag bundling and value concatenation**: clustered short
   flags (`-bS` == `-b -S`) and value-attached short flags (`-q20` == `-q 20`)
   are accepted, matching the getopt parsers of the upstream C tools

## POSIX short-flag bundling and value concatenation

Upstream bioinformatics tools parse their command lines with `getopt(3)`, which
accepts two short-flag idioms that Go's standard `flag` package does **not**
support out of the box:

- **Bundling** — several boolean short flags may share one dash:
  `samtools view -bS in.bam` is exactly `-b -S`. Any run of boolean short flags
  may be clustered: `-hb`, `-bSH`, etc.
- **Value concatenation** — a value-taking short flag may carry its value with
  no separating space: `-q20` is `-q 20`, `-@4` is `-@ 4`.
- **Mixed clusters** — bundling and concatenation combine left-to-right. The
  first value-taking flag in a cluster ends it and consumes the remainder as
  its value: `-bSq20` is `-b -S -q 20`. When that flag is the last character
  (`-bSq 20`), the value is taken from the following argument.

These idioms compose with the rest of the conventions: a long option
(`--min-mapq=20`), the `--` terminator, and a bare `-` (stdin/stdout) are all
passed through untouched. An unknown short character in a cluster is rejected
with a "flag provided but not defined" error (exit code 2), exactly as upstream
getopt rejects an unknown option.

To get this behaviour, parse the `*flag.FlagSet` through
[`cliflag.Parse`](../pkg/cliflag/posix.go) instead of calling `fs.Parse`
directly:

```go
// Drop-in replacement for fs.Parse(args). On error, print usage and exit 2.
if err := cliflag.Parse(fs, args); err != nil {
    fmt.Fprintln(os.Stderr, err)
    fmt.Fprint(os.Stderr, usage)
    return 2
}
```

`cliflag.Parse` introspects the FlagSet to distinguish boolean switches (whose
`flag.Value` reports `IsBoolFlag() == true`) from value-taking flags, expands
the clusters into Go-flag-canonical tokens, and then calls `fs.Parse`. It is
idempotent on already-canonical arguments, so wiring an existing tool through it
never changes the meaning of a command line that already worked. (Use
`cliflag.Normalize` if you only need the rewritten argument slice.)

## Option Naming Conventions

### Short Options (Single Letter)

Short options use a single dash followed by a single letter:

- `-i` for input files
- `-o` for output files
- `-q` for quality-related options
- `-l` for length-related options
- `-g` for GC content options
- `-t` for trimming options
- `-d` for duplicate removal options
- `-h` for help
- `-v` for version

### Long Options (Full Word)

Long options use a double dash followed by a descriptive word or hyphenated phrase:

- `--input` for input files
- `--output` for output files
- `--quality` for quality thresholds
- `--min-length` for minimum length
- `--max-length` for maximum length
- `--help` for help
- `--version` for version

### Paired Options

When tools support paired-end reads, use numbered suffixes:

- `-i1, --input1` for first read file
- `-i2, --input2` for second read file
- `-o1, --output1` for first output file
- `-o2, --output2` for second output file

## Standard Options

All tools should support these standard options:

- `-h, --help`: Display help information
- `-v, --version`: Display version information
- `-i, --input`: Primary input file (use `-` for stdin)
- `-o, --output`: Primary output file (use `-` for stdout if not specified)

## File Format Options

For tools that support multiple formats:

- `--fasta`: Input is FASTA format
- `--fastq`: Input is FASTQ format
- `--format`: Specify format explicitly (fasta, fastq, etc.)

## Common Filtering Options

Standard options across filtering tools:

- `-l, --min-length INT`: Minimum sequence length
- `-L, --max-length INT`: Maximum sequence length
- `-g, --min-gc FLOAT`: Minimum GC content percentage
- `-G, --max-gc FLOAT`: Maximum GC content percentage
- `-q, --min-quality FLOAT`: Minimum quality score
- `-Q, --max-quality FLOAT`: Maximum quality score
- `-n, --max-ns INT`: Maximum number of N bases
- `-N, --max-ns-percent FLOAT`: Maximum percentage of N bases

## Trimming Options

Standard trimming options:

- `--trim-left INT`: Trim N bases from 5' end
- `--trim-right INT`: Trim N bases from 3' end
- `--trim-qual-left INT`: Quality-based trimming from 5' end
- `--trim-qual-right INT`: Quality-based trimming from 3' end
- `--trim-n-left INT`: Trim poly-N from 5' end
- `--trim-n-right INT`: Trim poly-N from 3' end

## Duplicate Removal Options

- `-d, --derep MODE`: Duplicate removal mode (1=exact, 4=revcomp, 5=both)
- `--derep-min INT`: Minimum occurrences to keep

## Output Options

- `-o, --output FILE`: Main output file
- `--out-good FILE`: Passing sequences (alternative naming)
- `--out-bad FILE`: Rejected sequences
- `--stats FILE`: Statistics output file

## Help Text Format

Help text should follow this structure:

```
tool-name - Brief description

Usage:
  tool-name <command> [options]

Commands:
  command1    Description
  command2    Description

Options:
  -s, --short-name TYPE    Description
                           Additional description lines indented

Examples:
  tool-name command1 -i input.fastq -o output.fastq
  tool-name command2 --min-length 100 input.fasta
```

## Examples

### Good Example (PRINSEQ filter)

```bash
# Using short options
prinseq filter -i reads.fastq -o filtered.fastq -l 100 -q 20

# Using long options
prinseq filter --input reads.fastq --output filtered.fastq --min-length 100 --min-quality 20

# Mixed short and long options
prinseq filter -i reads.fastq --min-length 100 --trim-qual-left 20
```

### Good Example (paired-end)

```bash
# Short options
prinseq filter -i1 R1.fastq -i2 R2.fastq -o1 out_R1.fastq -o2 out_R2.fastq

# Long options
prinseq filter --input1 R1.fastq --input2 R2.fastq --output1 out_R1.fastq --output2 out_R2.fastq
```

## Implementation Guidelines

### Go Implementation

For Go tools using the standard `flag` package, register flags through the
shared [`pkg/cliflag`](../pkg/cliflag/cliflag.go) helpers. Each helper registers
the same destination under both the short and long name on one `*flag.FlagSet`,
so a single call gives you both forms. Either name may be empty to register only
one form (e.g. a long-only `--input-fmt-option`).

```go
fs := flag.NewFlagSet("mytool", flag.ContinueOnError)

cliflag.StringVar(fs, &input, "i", "input", "", "Input file (use '-' for stdin)")
cliflag.IntVar(fs, &minLen, "l", "min-length", 0, "Minimum sequence length")
cliflag.BoolVar(fs, &showHelp, "h", "help", false, "Show help")
```

Available helpers (all in `pkg/cliflag`):

| Helper                 | Backing type      | Typical use                         |
| ---------------------- | ----------------- | ----------------------------------- |
| `StringVar`            | `string`          | paths, format names                 |
| `IntVar`               | `int`             | counts, thresholds                  |
| `Int64Var`             | `int64`           | seeds, large counts                 |
| `Uint64Var`            | `uint64`          | non-negative large counts           |
| `Float64Var`           | `float64`         | quality/fraction thresholds         |
| `BoolVar`              | `bool`            | toggles, `--help`/`--version`       |
| `DurationVar`          | `time.Duration`   | timeouts                            |
| `Var`                  | `flag.Value`      | repeatable flags / custom parsing   |

Use `cliflag.Var` for repeatable flags (e.g. `-r/--region` that may appear
multiple times): both forms accumulate into the same destination. Only
hand-roll `fs.*Var` directly when a flag is intentionally single-form to match
an upstream tool (e.g. samtools' `-Q`/`-D` phase flags), and note why.

`cliflag.FormatUsage(short, long, valueType, description)` renders an aligned
`-s, --long TYPE   description` usage line.

### Validation

- Validate that required options are provided
- Provide clear error messages when options conflict
- Show relevant help text when errors occur

## Version History

- 2025-10-20: Initial version established for PRINSEQ tool
- 2026-05-13: Formalised POSIX-compliant CLI as a hard requirement for any
  tool to be considered "complete" (alongside 100% feature parity)
- 2026-06-10: Added `cliflag.Parse`/`cliflag.Normalize` for getopt-compatible
  POSIX short-flag bundling (`-bS`) and value concatenation (`-q20`); `samtools
  view` is the first tool wired through it, with repo-wide rollout to follow
- 2026-06-10: Rolled `cliflag.Parse` out across the htslib-family CLIs — all
  remaining `samtools` subcommands (sort, index, flagstat, depth, fastq,
  mpileup, idxstats, quickcheck, dict, cat, reheader, addreplacerg, fixmate,
  merge, coverage, split, markdup, stats, calmd, import, phase, targetcut,
  consensus), plus `tabix`, `bgzip`, and `htsfile`. Each subcommand now also
  registers the upstream getopt short flags it previously lacked (legacy /
  no-op compat flags) so bundled clusters that include them parse and behave
  like upstream. Repeatable `-a` (mpileup, depth) is now a count flag so the
  fused `-aa` ("all positions, all chromosomes") resolves through the same
  bundling path as upstream rather than a bespoke pre-pass.
- Future: May be extended as more tools are ported

## References

- POSIX Utility Conventions: <https://pubs.opengroup.org/onlinepubs/9699919799/basedefs/V1_chap12.html>
- GNU Coding Standards: <https://www.gnu.org/prep/standards/html_node/Command_002dLine-Interfaces.html>
