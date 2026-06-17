#!/usr/bin/env python3
"""Regenerate the FASTQ input fixtures for the fastp parity test.

The fixtures are checked in (so the parity test is reproducible without
needing Python at test time), but this script makes the contents
auditable. Run it from this directory:

    python3 generate.py

Every read is constructed deterministically (random.seed(42)).
"""
import random

random.seed(42)

ADAPTER = "AGATCGGAAGAGC"


def fmt(records):
    out = []
    for n, seq, qual in records:
        out.append(f"@{n}\n{seq}\n+\n{qual}\n")
    return "".join(out)


def write(path, records):
    with open(path, "w") as f:
        f.write(fmt(records))


def hi_q(n):
    return "".join(chr(33 + random.randint(30, 40)) for _ in range(n))


def med_q(n):
    return "".join(chr(33 + random.randint(20, 35)) for _ in range(n))


def lo_q(n):
    return "".join(chr(33 + random.randint(2, 10)) for _ in range(n))


def main():
    # 1) clean SE reads: 30 records, 100bp, high quality, random ACGT.
    clean = []
    for i in range(30):
        seq = "".join(random.choice("ACGT") for _ in range(100))
        clean.append((f"clean_{i}", seq, hi_q(100)))
    write("se_clean.fq", clean)

    # 2) SE reads with the canonical Illumina TruSeq adapter embedded mid-read.
    adapter = []
    for i in range(20):
        insert_len = random.randint(30, 60)
        insert = "".join(random.choice("ACGT") for _ in range(insert_len))
        rest = "".join(
            random.choice("ACGTN")
            for _ in range(100 - insert_len - len(ADAPTER))
        )
        seq = (insert + ADAPTER + rest)[:100]
        adapter.append((f"adp_{i}", seq, hi_q(100)))
    write("se_adapter.fq", adapter)

    # 3) Reads with N's scattered through; tests --n_base_limit boundary.
    ns = []
    for i in range(10):
        seq = "".join(
            random.choice("ACGTN") if random.random() < 0.3 else random.choice("ACGT")
            for _ in range(80)
        )
        ns.append((f"ns_{i}", seq, hi_q(80)))
    write("se_ns.fq", ns)

    # 4) Low-complexity mix: alternating AT-only and full-complexity reads.
    lowcomp = []
    for i in range(10):
        if i % 2 == 0:
            seq = ("AT" * 50)[:80]  # 100% alternation: complexity=1.0.
        else:
            seq = "".join(random.choice("ACGT") for _ in range(80))
        lowcomp.append((f"lc_{i}", seq, hi_q(80)))
    write("se_lowcomp.fq", lowcomp)

    # 5) Low-complexity homopolymer mix: pure homopolymers (complexity=0) plus
    # mixed reads, so the filter actually rejects something.
    homop = []
    for i in range(10):
        if i % 2 == 0:
            seq = ("A" * 80)  # 0% adjacent differences → complexity=0.
        else:
            seq = "".join(random.choice("ACGT") for _ in range(80))
        homop.append((f"hp_{i}", seq, hi_q(80)))
    write("se_homopolymer.fq", homop)

    # 6) Low-quality 3' tail reads (60 hi-Q + 40 lo-Q). For sliding-window
    # cut_right/cut_tail/cut_front parity.
    lq = []
    for i in range(15):
        seq = "".join(random.choice("ACGT") for _ in range(100))
        lq.append((f"lq_{i}", seq, hi_q(60) + lo_q(40)))
    write("se_lqtail.fq", lq)

    # 7) UMI reads: 8bp leading UMI + 72bp insert.
    umi = []
    for i in range(15):
        umi_seq = "".join(random.choice("ACGT") for _ in range(8))
        rest = "".join(random.choice("ACGT") for _ in range(72))
        umi.append((f"umi_{i}", umi_seq + rest, hi_q(80)))
    write("se_umi.fq", umi)

    # 8) Reads for Q20/Q30 filtering: a mix of fully-hi-Q and partially-lo-Q.
    qmix = []
    for i in range(20):
        seq = "".join(random.choice("ACGT") for _ in range(80))
        if i < 10:
            qual = hi_q(80)
        else:
            qual = med_q(40) + lo_q(40)
        qmix.append((f"qm_{i}", seq, qual))
    write("se_qmix.fq", qmix)

    # 9) PE: R1 and R2 are the reverse complement of each other where they
    # overlap (synthetic short-insert pairs).
    r1 = []
    r2 = []
    for i in range(40):
        seq = "".join(random.choice("ACGT") for _ in range(80))
        r1.append((f"pe_{i}", seq, hi_q(80)))
        comp = seq[::-1].translate(str.maketrans("ACGTN", "TGCAN"))
        r2.append((f"pe_{i}", comp, hi_q(80)))
    write("pe_r1.fq", r1)
    write("pe_r2.fq", r2)

    # 10) Reads explicitly suitable for duplication evaluation: 50 reads
    # where 20 are duplicates of one another.
    dup = []
    template = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"[:60]
    for i in range(50):
        if i < 30:
            seq = "".join(random.choice("ACGT") for _ in range(60))
        else:
            seq = template
        dup.append((f"d_{i}", seq, hi_q(60)))
    write("se_dup.fq", dup)

    # 11) Paired-end overlap-correction fixture: 30 pairs of 100bp reads with
    # an ~80bp overlap and one low-quality base error injected into R1 inside
    # the overlap, so `--correction` corrects exactly one base per pair.
    def revcomp(s):
        return s[::-1].translate(str.maketrans("ACGTN", "TGCAN"))

    random.seed(101)
    cr1, cr2 = [], []
    for i in range(30):
        insert = "".join(random.choice("ACGT") for _ in range(120))
        r1 = list(insert[:100])
        q1 = ["I"] * 100
        r2 = revcomp(insert[20:120])
        q2 = ["I"] * 100
        pos = random.randint(40, 90)
        r1[pos] = "A" if r1[pos] != "A" else "C"
        q1[pos] = "#"  # phred 2, low quality -> correctable
        cr1.append((f"cp_{i}", "".join(r1), "".join(q1)))
        cr2.append((f"cp_{i}", r2, "".join(q2)))
    write("corr_r1.fq", cr1)
    write("corr_r2.fq", cr2)

    # 12) Overrepresentation/split fixture: 6000 SE reads of 100bp with a 50bp
    # motif injected into one third, large enough to trigger upstream's
    # overrepresentation analysis and multi-file splitting.
    random.seed(202)
    motif = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"
    ora = []
    for i in range(6000):
        if i % 3 == 0:
            seq = motif + "".join(random.choice("ACGT") for _ in range(50))
        else:
            seq = "".join(random.choice("ACGT") for _ in range(100))
        ora.append((f"o_{i}", seq, "I" * 100))
    write("ora.fq", ora)

    # 13) SE adapter auto-detection fixture: 12000 reads of 100bp, above
    # upstream's evaluator gate (>= 10000 records, evaluator.cpp:344), so the
    # kmer/nucleotide-tree SE auto-detect path actually fires. Most reads are
    # read-through: a random insert followed by the canonical TruSeq adapter
    # AGATCGGAAGAGC and a short fixed adapter continuation, padded with random
    # bases. A minority are adapter-free so detection has to discriminate.
    #
    # Detection is a sampling-dependent heuristic, so the matching parity test
    # validates a SIMILARITY bound (per-read trimmed-length agreement +
    # base-identity) against the upstream binary, not strict byte-equality.
    random.seed(303)
    # The full Illumina TruSeq Read1 adapter; the leading 13bp
    # (AGATCGGAAGAGC) is the canonical known-adapter prefix fastp ships.
    full_adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
    detect = []
    for i in range(12000):
        if i % 10 == 0:
            # 10% adapter-free, full-length random reads.
            seq = "".join(random.choice("ACGT") for _ in range(100))
        else:
            insert_len = random.randint(20, 70)
            insert = "".join(random.choice("ACGT") for _ in range(insert_len))
            seq = (insert + full_adapter)
            if len(seq) < 100:
                seq = seq + "".join(
                    random.choice("ACGT") for _ in range(100 - len(seq))
                )
            seq = seq[:100]
        detect.append((f"det_{i}", seq, hi_q(100)))
    write("se_detect.fq", detect)


if __name__ == "__main__":
    main()
