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


//


// Fastq files
R1 = file(params.fastq1)
R2 = file(params.fastq2)


// --------------------------------------------------------------------------
// Step 0) Get data dir
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
// Step 1) Run FastQC on Raw Fastq Data
// --------------------------------------------------------------------------

process fastqc-00 {
  input:
      file FQ1 from R1
  input:
      file FQ2 from R2

  """
  fastqc $FQ1 $FQ2
  """
}

// --------------------------------------------------------------------------
// Step 2) Run Trimmomatic or other selected QC tool
// --------------------------------------------------------------------------

trim='trimmomatic'

process qc-fastq-01 {
  input:

  script:
  if( trim == 'trimmomatic')
  """
  """

  else if( trim == 'skewer')
  """
  """

  else if( trim == 'no')
  """
  """

  else
    error "No Trimming option selected: ${trim}"
}

// --------------------------------------------------------------------------
// Step 3) Run Alignment
// --------------------------------------------------------------------------

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
