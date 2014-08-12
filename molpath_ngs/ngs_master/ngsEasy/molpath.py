#!/usr/bin/env python

from ngsEasy.helpers import SGE

__doc__=='''
shell script based pipeline (adapted from R/bash scripts)
'''

def trimAdapters(p,sge):
    jobname = "novoalign_"+p.uniqueID()
    # determine script to run
    if p.RGPL().startswith('illumina'):
        if p.isPE():
            logging("Reads are PE ILLUMINA")
            runscript = "ngs_novoalign.illumina."+p.QUAL+".PE.sh"
        else:
            logging("Reads are SE ILLUMINA")
            runscript = "ngs_novoalign.illumina."+p.QUAL+".SE.sh"

    elif p.RGPL().startswith("iontorrent"):
        try:
            assert p.isPE()
        except AssertionError:
            raise Exception("IONTORRENT cannot be PE")
        else:
            logging("Reads are SE IONTORRENT")
            runscript = "ngs_novoalign.IONTORRENT."+p.QUAL+".SE.sh"
    else:
        raise Exception("mapping configuration????")

    # send job
    sge.qsub(exe=p.pipeline_dir+'/'+runscript+" "+p.fastq_prefix()+" "+p.uniqueID()+" "+p.wd(), name=jobname)

    # add dependency for next job
    sge.depends.append(jobname)

# read mapping (PE/SE, illumina/iontorrent)
def mapReads(p,sge):
    jobname = "novoalign_"+p.uniqueID()
    # determine script to run
    if p.RGPL().startswith('illumina'):
        if p.isPE():
            logging("Reads are PE ILLUMINA")
            runscript = "ngs_novoalign.illumina."+p.QUAL+".PE.sh"
        else:
            logging("Reads are SE ILLUMINA")
            runscript = "ngs_novoalign.illumina."+p.QUAL+".SE.sh"

    elif p.RGPL().startswith("iontorrent"):
        try:
            assert p.isPE()
        except AssertionError:
            raise Exception("IONTORRENT cannot be PE")
        else:
            logging("Reads are SE IONTORRENT")
            runscript = "ngs_novoalign.IONTORRENT."+p.QUAL+".SE.sh"
    else:
        raise Exception("mapping configuration????")

    # send job
    sge.qsub(exe=p.pipeline_dir+'/'+runscript+" "+p.fastq_prefix()+" "+p.uniqueID()+" "+p.wd(), name=jobname)

    # add dependency for next job
    sge.depends.append(jobname)

    return

def sam2bam(p,sge):
    jobname = 'sam2bam.'+p.uniqueID()
    runscript = 'ngs_sam2bam.sh'
    sge.qsub(exe=p.pipeline_dir+'/'+runscript+" "+p.uniqueID()+" "+p.wd(), name=jobname)
    sge.depends.append(jobname)
    return

# sort reads
def sortSam(p, sge):
    sge.qsub(exe=p.pipeline_dir+'/'+runscript+" "+p.uniqueID()+" "+p.wd()+" "+p.tmp(), name=jobname)
    sge.depends.append(jobname)
    return

${ngs_pipeline}/ngs_sam2bam.sh ${sample_name} ${sample_dir};

${ngs_pipeline}/ngs_SortSam.sh ${sample_name} ${sample_dir} ${sample_temp};
