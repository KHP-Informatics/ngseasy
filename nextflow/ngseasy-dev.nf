#!/usr/bin/env nextflow

/*
 * Main NGSeasy pipeline script
 *
 * @authors
 * Stephen Newhouse <stephen.j.newhouse@gmail.com>
 *
 */

 // nextflow run ngseasy.nf --studyId --sampleID --projectDir --

// define run parameters

params.studyId =
params.sampleId =
params.ProjectDir =
params.fastq1 =
params.fastq2 =
params.refgenome =
params.fastqc = 
params.trimmomatic =
params.aligner =
params.indelRealigner =
params.baserecal =
params.variantCaller =
params.svCaller =




process run-fastqc-00 {

  """
  fastqc
  """
}

process run-trimmomatic-01 {
}

process run-alignment {
}

process ngs-realign {
}

process ngs-bsqr {
}

process ngs-varcall {
}

process ngs-svcall {
}
