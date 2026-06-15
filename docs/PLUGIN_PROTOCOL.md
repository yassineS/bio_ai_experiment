# bcftools plugin protocol

This document specifies the contract between the Go `bcftools` host and a
*plugin*. It applies to the `bcftools plugin <name>` subcommand and the
`bcftools +<name>` shorthand.

## Why this differs from upstream

Upstream bcftools loads plugins as native shared objects (`.so`) via
`dlopen`, and a plugin must export a fixed C ABI (`init`, `process`,
`version`, ...). That model ties every plugin to a C toolchain, a matching
bcftools build, and the host's address space.

This port deliberately uses a **subprocess protocol** instead: a plugin is an
ordinary executable that the host runs as a child process, streaming variant
data over standard streams. A plugin can therefore be written in *any*
language — Go, Python, Rust, awk, a shell script — and needs no headers, no
ABI version check, and no rebuild when the host changes. The trade-off is the
serialization cost of piping VCF text; for the experiment's goals (let users
write their own plugins) that is the right call.

In addition to this subprocess protocol, the host now bundles **native
(pure-Go) reimplementations of all 41 upstream in-tree plugins**, dispatched by
`+<name>` ahead of the subprocess lookup. The subprocess protocol described in
this document is the **fallback** for any name not in the native registry — so
user-supplied executables in any language keep working unchanged. Native
plugins implement the `NativePlugin` interface in
`tools/bcftools/pkg/bcftools/native_plugin.go` (an in-process `Init`/`Process`/
`Destroy` contract mirroring upstream's plugin lifecycle, with optional
parallel/buffered/run-style/multi-output extensions), not this stdin/stdout
contract; a single reference subprocess plugin still ships as an example and
test fixture (`tools/bcftools/plugins/example`).

## Discovery

The host resolves a plugin name to an executable using the
`BCFTOOLS_PLUGINS` environment variable, a colon-separated (`os.PathListSeparator`)
list of directories — the same variable upstream uses.

- `bcftools +<name>` (or `bcftools plugin <name>`) looks for a file named
  **exactly `<name>`** in each `BCFTOOLS_PLUGINS` directory, in order, and
  uses the first one that is a regular file with an executable bit set. The
  first directory wins.
- A leading `+` on the requested name is stripped before lookup.
- If `<name>` contains a path separator it is treated as an explicit path and
  the search directories are bypassed.
- If `BCFTOOLS_PLUGINS` is unset or empty, no plugins are discoverable.

`bcftools plugin -l` (or `--list-plugins`) lists the discoverable plugin
names, one per line, de-duplicated and sorted. `bcftools plugin -lv` lists
them verbosely: the absolute path of each plugin, optionally followed by a
tab and the plugin's one-line `--about` description (see below).

## Execution model

`bcftools +<name> [plugin-args...] [-- <input> [regions...]]` runs as follows:

1. **Input normalisation (host side).** The host reads its input
   (`<input>`, or stdin when omitted/`-`) — VCF, VCF.gz, or BCF — applies any
   region/`-R` selection, and serialises the result to **uncompressed VCF
   text**.
2. **Plugin invocation.** The host executes the resolved plugin with
   `[plugin-args...]` as its `argv` (everything between the plugin name and a
   literal `--`). It pipes the VCF text from step 1 to the plugin's **stdin**
   and captures the plugin's **stdout**.
3. **Output formatting (host side).** The host parses the plugin's stdout as
   VCF and re-emits it in the container requested by `-O/--output-type`
   (`v` = VCF, `z` = VCF.gz; `b`/`u` once the BCF writer lands), writing to
   `-o/--output` or stdout.

So, from the plugin author's point of view:

> **A plugin is a filter that reads a VCF on stdin and writes a VCF on stdout.**

### Streams

| Stream | Direction | Format |
| ------ | --------- | ------ |
| stdin  | host -> plugin | uncompressed VCF text (header + records) |
| stdout | plugin -> host | uncompressed VCF text (header + records) |
| stderr | plugin -> host | free-form; forwarded verbatim to the host's stderr |

The plugin **must** emit a complete, valid VCF on stdout — at minimum the
`##fileformat` line, the `#CHROM` header line, and zero or more data records.
A plugin that filters records still emits the (possibly edited) header.

### argv

Everything on the command line after the plugin name and before a literal
`--` is passed verbatim as the plugin's `argv` (after `argv[0]`). Host
options (`-o`, `-O`, `-r`, `-R`, `-l`, ...) must appear *before* the plugin
name. The `--` separator, when present, ends the plugin arguments; the tokens
after it are the host input file and optional regions.

```
bcftools plugin myplugin --tag AF --threshold 0.01 -- input.vcf.gz chr1:1-1000
                 \_____/ \______________________/    \____________________/
                  name        plugin argv               host input + region
```

### Exit codes

- A plugin that exits **0** is considered successful; its stdout is used.
- A plugin that exits **non-zero** is an error: the host reports a
  `plugin "<name>" failed` error (exit status 1) and surfaces the plugin's
  captured stderr.
- If the plugin cannot be found or is not executable, the host fails before
  running anything.

### The optional `--about` flag

`bcftools plugin -lv` probes each plugin by running it with a single
`--about` argument. A plugin that supports it should print a one-line
human-readable description to stdout and exit 0. This is **optional**: a
plugin that does not recognise `--about` (and instead, say, tries to read a
VCF from stdin) is still listed — it simply shows no description. Plugin
authors who want a nice `-lv` entry should special-case `--about` early,
before touching stdin.

## Writing a plugin

Any of these is a valid plugin:

```sh
#!/bin/sh
# passthru: identity filter
exec cat
```

```python
#!/usr/bin/env python3
import sys
if "--about" in sys.argv[1:]:
    print("drop-id: blank out the ID column"); sys.exit(0)
for line in sys.stdin:
    if line.startswith("#"):
        sys.stdout.write(line)
    else:
        f = line.rstrip("\n").split("\t")
        f[2] = "."            # ID is column 3
        sys.stdout.write("\t".join(f) + "\n")
```

See `tools/bcftools/plugins/example/main.go` for a documented Go example that
resets every record's FILTER to `PASS` and reports a record count on stderr.

## Limitations

- Variant data crosses the process boundary as VCF text, so BCF-only binary
  fidelity is not preserved through a plugin; the host re-parses the plugin's
  VCF output before applying `-O`.
- The plugin sees the *normalised* VCF (post region selection), not the
  original file bytes.
- There is no streaming back-pressure contract beyond ordinary OS pipe
  buffering; very large inputs are buffered in memory by the host in v1.
