#########################################################################
# -- Author: Amos Folarin                                               #
# -- Organisation: KCL/SLaM                                             #
# -- Email: amosfolarin@gmail.com                                       #
#########################################################################


#------------------------------------------------------------------------
# This is the main calling script for scheduling sge jobs for the cnv pipeline
# INPUT A) FASTQs B) BAMs depending on what we get from the sequencing centre 
# The current sequence of steps will be something like:
# A.1) initial qc: macs, fastqc, fastx_clipper, fastq_quality_trimmer, fastx_artifacts_filter, fastq_quality_filter  etc.
# A.2) alignment: novoalign
# B.1) BAM: remove duplicates, base quality score recalibration, INDEL realign
# B.2) GATK cleaning
# C.1) CNV calling: GenomeStripe, Pindel
# )
# )
# )
# )
#------------------------------------------------------------------------


#------------------------------------------------------------------------
# 
#------------------------------------------------------------------------
