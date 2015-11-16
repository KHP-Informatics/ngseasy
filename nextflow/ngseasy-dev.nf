#!/usr/bin/env nextflow

/*
 * Main NGSeasy pipeline script
 *
 * @authors
 * Stephen Newhouse <stephen.j.newhouse@gmail.com>
 *
 */

// define run parameters

params.studyId =
params.sampleId =
params.ngsProjectDir =
params.fq1 =
params.fq2 =
params.refgenome =
params.fastqc =
params.trimmomatic =
params.aligner =
params.indelRealigner =
params.baserecal =
params.variantCaller =
params.svCaller =

process fastqQc {
}

process qctrimm {
}

process align {
}

process realign {
}

process bqqr {
}

process varcall {
}

process svcall {
}
