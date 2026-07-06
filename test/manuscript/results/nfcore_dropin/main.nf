#!/usr/bin/env nextflow
// Minimal DSL2 wrapper that runs the REAL, UNMODIFIED nf-core samtools/flagstat
// module (vendored verbatim under modules/nf-core/samtools/flagstat/main.nf).
// The point of the demo is to prove our Go `samtools` is a drop-in for the
// module: this wrapper is the only thing we add — the module file itself is
// byte-for-byte the upstream nf-core one. Which `samtools` runs is decided
// solely by PATH (see run_dropin.sh), and containers/conda are disabled in
// nextflow.config so the local executor resolves `samtools` from PATH.
nextflow.enable.dsl = 2

include { SAMTOOLS_FLAGSTAT } from './modules/nf-core/samtools/flagstat/main.nf'

workflow {
    ch_input = Channel.of(
        [ [ id: 'test' ], file("${projectDir}/data/test.bam"), file("${projectDir}/data/test.bam.bai") ]
    )
    SAMTOOLS_FLAGSTAT(ch_input)
    SAMTOOLS_FLAGSTAT.out.flagstat.view { meta, f -> "flagstat -> ${f}" }
}
