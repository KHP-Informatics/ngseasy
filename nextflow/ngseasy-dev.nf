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

process ngs-fastqc {
  container 'compbio/ngseasy-fastqc:v1.0-r002'

  """
  hshshs
  """
}

process ngs-qctrimm {
}

process ngs-align {
}

process ngs-realign {
}

process ngs-bsqr {
}

process ngs-varcall {
}

process ngs-svcall {
}
