#!/usr/bin/env python

from ruffus import *

__doc__=='''
ruffus based pipeline steps
'''

def NOTREADY():
    raise Exception('Not implemented yet')


def trim_adapters():
    qsub -o ${SGE_OUT} -e ${SGE_OUT} -q ${queue_name} -N novoalign.${sample_name} -l h_vmem=${novo_mem}G -pe multi_thread ${novo_cpu} -M ${email_contact} -m beas ${ngs_pipeline}/ngs_novoalign.illumina.${qual_type}.PE.sh ${fastq_prefix} ${sample_name} ${sample_dir};

def map_sequence
    NOTREADY()

def