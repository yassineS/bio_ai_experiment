# 03 — Related work & annotated bibliography

Positioning + a verification-graded reference list, synthesized from the research briefs.

**Verification grades:** `[V]` corroborated across ≥2 independent sources;
`[S]` snippet-sourced (exact number not read from primary PDF — **verify before print**);
`[PR]` vendor/CEO/press claim (not peer-reviewed); `[26]` 2026-dated arXiv id to re-confirm.
Many primary PDFs were egress-blocked during research; treat every specific percentage as
**[verify-PDF]** until checked.

---

## 1. Enterprise LLM legacy modernization (the closest "it works at scale" prior art)

- **Nikolov et al., "How is Google using AI for internal code migrations?"** arXiv:2501.06972,
  ICSE-SEIP 2025. ~80% of landed changes AI-authored; ~50% time reduction; int32→int64 across
  500M-LOC Ads codebase; JUnit3→4 (5,359 files / ~149k LOC / 3 mo). `[V]`/`[S]`
- **Ziftci et al., "Migrating Code At Scale With LLMs At Google."** arXiv:2504.09691, FSE 2025
  (Industry). 39 migrations/12 mo; 595 changes / 93,574 edits; 74.45% of changes & 69.46% of
  edits LLM-generated. **Best peer-reviewed enterprise data point.** `[V]`/`[S]`
- **Amazon Q Code Transformation (Jassy, 2024).** ~30k apps Java 8/11→17; "4,500 developer-years",
  "$260M", 79% shipped unchanged. **`[PR]` — unaudited; one customer reported ~36% acceleration
  (The Register).** Do **not** cite "260,000 dev-years/hours" (a conflation).
- **IBM watsonx Code Assistant for Z** (COBOL→Java). `[PR]`, qualitative only; the rigorous
  adjacent item is the eval-methodology preprint *Quality Evaluation of COBOL to Java Code
  Transformation* (arXiv:2507.23356).
- **XMainframe** (arXiv:2408.04660), **MainframeBench**, **COBOLEval** — COBOL LLM benchmarks. `[V]`

*Positioning:* all are **same-or-cross-language migration of proprietary code, test-suite-gated,
counterfactual self-estimated.** None is whole-tool byte-exact parity on public tools.

## 2. The LLM-coding evaluation crisis (why our oracle matters)

**Benchmarks:** SWE-bench (Jimenez, arXiv:2310.06770, ICLR'24; 2,294 inst./12 Python repos;
Claude 2 ~1.96% BM25) `[V]`; SWE-bench Verified (OpenAI 2024, 500 inst.) `[V]`; Multimodal
(arXiv:2410.03859, ~617 inst./17 JS repos) `[V]`; Multilingual (300/42/9 langs); Pro (Scale,
arXiv:2509.16941, 1,865 inst., contamination-resistant) `[V]`; SWE-agent (arXiv:2405.15793,
NeurIPS'24, the ACI) `[V]`; OpenHands (arXiv:2407.16741, ICLR'25) `[V]`; Aider Polyglot;
Devin (Cognition 2024, 13.86% + "Debunking Devin") `[V]`/contested.

**Contamination:** SWE-Bench Illusion (arXiv:2506.12286; buggy-file-from-issue 76% in-bench vs
53% out; >94% predate cutoffs) `[V]`/`[S]`; *Does SWE-bench Test Ability or Memory?*
(arXiv:2512.10218; 3×/6×) `[S][26]`; SWE-rebench (arXiv:2505.20411); EvalPlus (arXiv:2305.01210,
NeurIPS'23; +80×/35× tests, −19–29% pass@k) `[V]`/`[S]`; Riddell et al. (ACL'24; 50.8/63.4%
Stack overlap); Roberts et al. (ICLR'24, cutoff natural experiment); LiveCodeBench (arXiv:2403.07974).

**Weak-oracle / "passes tests ≠ correct":** SWE-Bench+ (arXiv:2410.06992; 32.67% leakage, 31.08%
weak tests, 12.47%→3.97%) `[S]`; *Are "Solved Issues" Really Solved?* / PatchDiff (arXiv:2503.15223,
ICSE'26; ~29.6% suspicious) `[S]`; UTBoost (arXiv:2506.09289; re-labels 24–41%) `[S]`; OpenAI
"Why we no longer evaluate SWE-bench Verified" (2026) `[V]`.

**Reward hacking / spec gaming (what the external-binary oracle prevents):** Anthropic
*Natural Emergent Misalignment from Reward Hacking* (arXiv:2511.18397; `AlwaysEqual`/`sys.exit(0)`/
`conftest.py`, ~12% sabotage) `[V]`/`[S]`; OpenAI CoT monitoring (arXiv:2503.11926;
`exit(0)`/`SkipTest`, obfuscation risk) `[V]`; METR (2025; scorer-clock/answer-exfiltration,
43× RE-Bench) `[V]`/`[S]`; ImpossibleBench (arXiv:2510.20270; GPT-5 76% cheat) `[S]`; EvilGenie
(arXiv:2511.21654) `[26]`; Bondarenko (arXiv:2502.13295; chess board-file edit) `[V]`; Krakovna
spec-gaming list (2018, GenProg empty-list, file-deletion) `[V]`.

**Package hallucination / supply chain:** Spracklen, *We Have a Package for You!* (USENIX'25,
arXiv:2406.10279; 19.7% hallucinated, 21.7% open-source, 205,474 names, 43% reproducible →
slopsquatting) `[V]`/`[S]`.

## 3. Differential testing, the oracle problem, APR overfitting (our theoretical backbone)

**Oracle problem & differential lineage:** Barr et al., *The Oracle Problem* (IEEE TSE 2015;
pseudo-oracle = N-version/reference-implementation oracle) `[V]`; McKeeman 1998 (coins
differential testing); Csmith (PLDI'11, >325 compiler bugs) `[V]`; EMI (PLDI'14, 147 bugs) `[V]`;
NEZHA (S&P'17, domain-independent) `[V]`.

**Reference-implementation-as-oracle (closest method match):** SemGraft (Mechtaev, ICSE'18;
BusyBox vs Coreutils) `[V]`; JEST (arXiv:2102.07498, N+1-version JS engines); Mokav
(arXiv:2406.10375, JSS'25; difference-exposing tests); **RustAssure** (arXiv:2510.07604;
original C as oracle for LLM-transpiled Rust — *the* methodological cousin) `[V]`.

**APR overfitting ("plausible ≠ correct"):** Z. Qi et al. (ISSTA'15; Kali, GenProg correct on
**2/105**, plausible-but-incorrect dominates) `[V]`; Long & Rinard (ICSE'16; plausible patches
orders of magnitude denser than correct) `[V]`; Smith et al. (FSE'15; coins "overfitting",
held-out IntroClass) `[V]`; Monperrus bibliography (CSUR'18) `[V]`. *Corrections vs. common
mis-citation:* Kali = **Zichao** Qi ISSTA'15 (not Yuhua Qi ICSE'14); the "Critical Review" is
**Monperrus solo, ICSE'14** (arXiv:1408.2103); ODS is **IEEE TSE 2022**.

*Positioning:* our method is the APR/oracle literature's prescribed answer (reference-as-oracle,
N-version differential testing) **scaled from function/patch to whole interdependent tools**,
with the original binary as oracle — exactly what RustAssure does for fragments, but for entire
shipped CLIs validated byte-for-byte.

## 4. AI rewrites of real tools — the novelty verdict ★

**No peer-reviewed byte-exact whole-tool parity result exists (mid-2026).** Closest cases:

- **wedeo** (Claude-built FFmpeg-in-Rust): bit-exact CI on *individual codecs* (H.264 79/79),
  but partial, slower, no whole tool. `[V]` (GitHub/HN)
- **RustQC** (Seqera, Claude-built): 15 QC tools→1 binary, **byte-identical featureCounts only**,
  self-reported, flagged experimental. `[V]`/`[PR]`
- **MirrorCode/Gotree** (METR+Epoch 2026): black-box reconstruction of 16.9k-LOC Gotree, graded
  by hidden tests — **not a documented byte-exact criterion.** `[V]`/`[26]`
- **Meta ProgramBench** (arXiv:2605.03546, 2026 preprint): whole-program reconstruction vs a
  reference binary via differential fuzzing over ~248k tests / 200 tasks (jq, ripgrep, FFmpeg,
  SQLite, …) — **no model fully solved any task.** `[S][26]`
- **Commit0** (arXiv:2412.01769, ICLR'25): from-scratch libraries graded by *unit tests*; no
  agent fully reproduces full libraries. `[V]`
- **uutils coreutils** (human Rust port): mature, Ubuntu-default, yet **50.6% of fuzzed `test`
  cases diverge** from GNU (Li et al., NDSS'25). `[V]` — the human-effort/difficulty baseline.
- **Heng Li, "The AI Rewrite Dilemma"** (2026): the samtools/bwa/minimap2 author's skeptical
  commentary; references RustQC; names validation+maintenance as the unsolved problem. *Ideal
  skeptical citation from inside the community.* `[V]`/`[26]`

*This is the wedge:* public, ubiquitous, interdependent scientific tools, byte-exact vs the
original binary, memory-safe/cgo-free, peer-reviewable and re-verifiable — the empty quadrant.

## 5. Bioinformatics domain: motivation, truth sets, validation precedent, integration

**Format specs / conformance corpora (the oracle substrate):** `samtools/hts-specs` (SAM v1.6,
CRAM v2.1/v3/codecs, VCF v4.1–4.5, BCF, BED, CSI/Tabix); GFF3 v1.26 (Sequence Ontology); htslib
`test/` fixtures; **htscodecs-corpus** (rANS/fqzcomp vectors); **OSS-Fuzz** `hts_open_fuzzer`
(never-crash standard; real `vcf_parse_format` crash). Originating papers: Li 2009 (SAM,
btp352); Danecek 2011 (VCF, btr330); Hsi-Yang Fritz 2011 (CRAM, gr.114819.110); Quinlan & Hall
2010 (BEDTools, btq033); **Danecek 2021 (HTSlib/Twelve years of SAMtools and BCFtools,
GigaScience giab008/giab007).** `[V]`

**Re-implementation precedent (for D1 framing):** **Sentieon = GATK
reimplementation, >99.99% concordant after removing GATK downsampling** (Kendig 2019, Front.
Genet.; Freed 2017) — the template of validating a re-implementation against the original tool
it replaces. `[V]`/`[S]`

**Pipeline integration / "drop-in" substrate (for C4):** Nextflow (Di Tommaso 2017, nbt.3820),
nf-core (Ewels 2020, s41587-020-0439-x), Snakemake (Köster 2012, bts480; Mölder 2021), Galaxy
(2024 update, gkae410), GATK Best Practices (McKenna 2010; DePristo 2011), Bioconda (Grüning
2018, s41592-018-0046-7), BioContainers (da Veiga Leprevost 2017, btx192). *Every orchestrator
addresses tools by literal CLI string + paths → flag+output identity = transparent swap.* `[V]`

## 6. Evaluation-norms references (for C3/C5 rigor)

Cost-aware eval: Kapoor et al., *AI Agents That Matter* (arXiv:2407.01502); HAL leaderboard
(arXiv:2510.11977). Nondeterminism: Thinking Machines *Defeating Nondeterminism in LLM Inference*
(2025, blog); *On Randomness in Agentic Evals* (arXiv:2602.07150) `[26]`; pass@k vs pass^k
(τ-bench, arXiv:2406.12045). Autonomy/HITL: Levels of Autonomy (Feng, arXiv:2506.12469); Human
Agency Scale H1–H5 (Stanford, arXiv:2506.06576); Levels of AGI autonomy (Morris, arXiv:2311.02462,
ICML'24). Reproducibility: ACM Artifact Review & Badging. Metamorphic testing in bioinformatics
(BMC Bioinformatics 1471-2105-10-24; MET'22). `[V]`/`[S]`

---

## Domain "what's hard" specifics (reviewer-probe list — see also `02 §B`)

Ranked by likelihood of silently breaking byte-exactness (for the failure-modes section):

1. **Float `%g` formatting + QUAL/PL last-ULP rounding** — Go `strconv` vs C `printf` differ in
   shortest-float digits and `inf`/`nan` spelling; glibc vs Go `math` differ at the last ULP.
2. **CRAM reference-MD5 (M5) handling** — wrong ref → *silent* base corruption (worst class);
   `REF_PATH`/`REF_CACHE`/ENA fallback must match.
3. **BGZF block boundaries → virtual offsets** — `.bai/.csi/.tbi` byte-exactness needs identical
   block-packing (why the deflate backend is pinned).
4. **Sort comparator + stability** — `strnum_cmp` natural ordering, `*`/unmapped placement,
   tie-break (in)stability (define byte-exactness at `-@1`).
5. **`bcftools norm` INFO/FORMAT Number=A/R/G re-indexing** on multiallelic split (Tan 2015).
6. **MD/NM edge cases** (`=`/`X`/`N`/ambiguity codes in calmd).

Plus structural genomics gotchas: 0-based BED vs 1-based VCF/SAM/GFF; chromosome naming/order
(`chr1` vs `1`); `LC_ALL=C` collation. "Byte-exact vs upstream" is only well-defined **against a
pinned upstream version** (the repo vendors them — cite the pinned commit).

---

## Citation hygiene before submission

- Re-pull every `[S]`/`[26]` figure from the primary PDF (research ran under egress blocks).
- Use the **vendored** specs/manuals as version-pinned citations (e.g.
  `reference_code/bcftools/doc/bcftools.txt` for `norm` semantics) rather than live URLs.
- Quarantine `[PR]` vendor numbers explicitly as industry claims without independent audit.
- Carry the three APR mis-citation corrections (Kali=Z.Qi ISSTA'15; Critical Review=Monperrus
  ICSE'14; ODS=TSE'22).
