#!/usr/bin/env nextflow

params.fq1 =
params.fq2 =  
params.refgenome =

process align {

	input:
	file

	output:
	file raw_bam

	"""
	docker run -it compbio/bwa > raw_bam 
	"""

}

