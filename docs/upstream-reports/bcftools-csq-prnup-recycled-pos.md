# Upstream bug report (ready to file): bcftools csq `@<pos>` marker prints a recycled record's position

This is a prepared bug report for the upstream **samtools/bcftools** project.
It is staged here (rather than filed automatically) because this session's
GitHub access is scoped to `yassineS/bio_ai_experiment`; please paste the body
below into a new issue at <https://github.com/samtools/bcftools/issues/new>.

Cross-reference: `docs/UPSTREAM_BUGS.md#bcftools-csq-prnup-recycled-pos`.

---

## Title

`bcftools csq`: compound `@<pos>` marker can print an unrelated (recycled) record's position

## Version

- bcftools 1.23.1 (`csq --force -f ref.fa -g annot.gff3 -p a in.vcf`), htslib 1.23.1.
- Reproduced from current `csq.c` / `vbuf` source; the relevant code is unchanged on the development branch.

## Summary

When a variant is a *silent member* of a haplotype-combined (compound)
consequence, `csq` prints a marker `@<pos>` whose `<pos>` is meant to be the
position of the record that carries the full consequence string (the haplotype
anchor). In some configurations `<pos>` is instead the position of a completely
unrelated, later record. The consequence is correctly attributed to the right
transcript/gene; only the printed anchor position is wrong.

## Root cause

The marker stores a **`bcf1_t*` pointer** to the anchor record and dereferences
`->pos` only later, at output time:

- `hap_add_csq()` (`csq.c`): for each silent member it sets
  `tmp_csq->type.ref = hap->stack[ref_node].node->rec;` — a pointer to the
  anchor record's `bcf1_t`.
- `kput_vcsq()` (`csq.c`): emits `@`, then `kputw(csq->ref->pos+1, str)` —
  reading `ref->pos` at print time.

Those `bcf1_t` records live in a **recycled ring buffer**. `vbuf_push()` does:

```c
if ( !vrec->line ) vrec->line = bcf_init1();
SWAP(bcf1_t*, (*rec_ptr), vrec->line);
```

so when a `vbuf` slot is reused for a later incoming record, the stored
`bcf1_t` (still pointed to by a haplotype node's `->rec`) is swapped back out
and refilled with the new record's data. If the anchor record's slot is
flushed and recycled **before** the silent member's `vbuf` is flushed,
`ref->pos` reads the recycled record's position.

This happens naturally when a variant is shared between the compound
haplotypes of **two overlapping genes**: one gene's anchor can be output (and
its slot recycled) while the silent member is still buffered, held back by the
*other* gene's still-active transcript.

## Observed vs. expected

On a fixture with two overlapping genes — `gene00024` (`+`, frameshift anchor
at chr1:30420) and `gene00025` (`-`, frameshift anchor at chr1:31623) — the
variants chr1:30690 / 30722 / 31112 are silent members of **both** frameshift
haplotypes. Output:

```
chr1  30690  ...  BCSQ=@33495,@31623   # observed
chr1  30690  ...  BCSQ=@30420,@31623   # expected
```

`@31623` (the gene00025 anchor) is correct. The first marker should be
`@30420` (the gene00024 anchor) but prints `@33495`, which is an unrelated
`intron|gene00025` record at chr1:33495 that reused the recycled slot.

Instrumenting `hap_add_csq()` confirms the marker is correctly attributed to
gene00024 with `ref->pos+1 == 30420` **at staging time**; instrumenting the
`@`-print in `kput_vcsq()` shows the same pointer reads `33495` **at output
time** — i.e. the pointer is stable but the `bcf1_t` behind it was recycled.

```
PRNUP var=30690 gene=gene00024 strand=1 refpos=30420 ibeg=1 iend=5 i=3 refnode=1   # staging
OUT   @33495 gene=gene00024                                                        # output
```

## Suggested fix

Resolve the anchor position by value when the marker is staged, rather than
holding a `bcf1_t*` and dereferencing it after the buffer may have recycled it.
For example, store `int32_t ref_pos` (and, if needed, the rid) in the `vcsq_t`
alongside / instead of `bcf1_t *ref`, set it in `hap_add_csq()`, and print that
in `kput_vcsq()`. Alternatively, ensure a silent member's referenced anchor is
never flushed before the member itself (keep the anchor pinned until all
records that reference it are emitted).

## Notes

Found while building an independent Go re-implementation of `csq` and
validating it byte-for-byte against the upstream binary. Our port resolves the
anchor position by value and therefore emits the expected `@30420`; everything
else in the haplotype-aware output (INFO/BCSQ consequence set and order, the
FORMAT/BCSQ sample bitmask, intron gene selection, splice-region ordering) is
byte-identical to upstream on this fixture.
