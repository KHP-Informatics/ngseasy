#!/usr/bin/env python

#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -m beas
#$ -pe multi_thread 1
#$ -l h_vmem=1G
#$ -p -0.99999999999999999999999999999999999999999999999999999999999999999
#$ -j y
#------------------------------------------------------------------------#

__doc__ = '''
#############################################################################################
# -- Authors: Stephen Newhouse, Amos Folarin, Aditi Gulati, David Brawand                                  #
# -- Organisation: KCL/SLaM/NHS                                                             #
# -- Email: stephen.j.newhouse@gmail.com, amosfolarin@gmail.com,aditi. gulati@nhs.net, dbrawand@nhs.net         #
# -- Verion: 1.3                                                                            #
# -- Date: 11/09/2013                                                                       #
# -- DESC: NGS pipeline to perform SE/PE Alignments & GATK cleaning                         #
#############################################################################################
'''
__author__ = "Stephen Newhouse, Amos Folarin, Aditi Gulati, David Brawand"
__copyright__ = ""
__credits__ = []
__license__ = "LGPL"
__version__ = "1.4"
__maintainer__ = "David Brawand"
__email__ = "stephen.j.newhouse@gmail.com, amosfolarin@gmail.com, aditi.gulati@nhs.net, dbrawand@nhs.net"
__status__ = "Development"  # ["Prototype", "Development",  "Production"]


import sys
import os
import re
import json
import logging

from optparse import OptionParser

from ngsEasy.Patient import Patients
from ngsEasy.helpers import syscall, sgecall
#import ngsEasy.pipeline as mp  # adaptation as a ruffus pipeline (incompatible with the docker containerisation)
import ngsEasy.molpath as mp  # the old school bash script based pipeline (all submitted at once)


# setup pipeline environment (this has to be streamlined)
def setupPipeline(config):
    # set own path
    os.environ["ngs_pipeline"] = os.path.abspath(os.path.dirname(sys.argv[0]))

    # tmp alignment dir
    os.environ["ngstmp"] = config['paths']['ngstmp']

    # SGE
    os.environ["queue_name"] = config['gridEngine']['queue']

    # tools
    os.environ["ngs_picard"] = config['software']['picard']
    os.environ["ngs_gatk"] = config['software']['gatk']
    os.environ["ngs_novo"] = config['software']['novo']
    os.environ["ngs_samtools"] = config['software']['samtools']
    os.environ["annovar"] = config['software']['annovar']
    os.environ["java_v1_7"] = config['software']['java_v1_7']

    # Resources
    # genome/reference
    os.environ["reference_genome_novoindex"] = config['reference']['genome']['novoindex']
    os.environ["reference_genome_seq"] = config['reference']['genome']['sequence']
    # known indel covariates for base recalibration
    os.environ["b37_Mills_Devine_2hit_indels"] = config['reference']['indels'][0]
    os.environ["b37_1000G_indels"] = config['reference']['indels'][1]
    # known SNP covariates for base recalibration
    os.environ["b37_1000G_omni2_5"] = config['reference']['snps'][0]
    os.environ["b37_dbsnp"] = config['reference']['snps'][1]
    os.environ["b37_hapmap_3_3"] = config['reference']['snps'][2]
    os.environ["b37_1000G_snps"] = config['reference']['snps'][3]
    os.environ["annovar_humandb"] = config['reference']['annotation']['annovar_humandb']

    # mem and cpu vars (TO BE RENAMED AND CLEANED UP)
    os.environ["novo_cpu"] = str(config['resources']['novo']['cpu'])
    os.environ["novo_mem"] = str(config['resources']['novo']['mem'])
    os.environ["sge_h_vmem"] = str(config['resources']['picard']['sge_h_vmem'])
    os.environ["java_mem"] = str(config['resources']['picard']['java_mem'])
    os.environ["gatk_h_vmem"] = str(config['resources']['gatk']['gatk_h_vmem'])
    os.environ["gatk_java_mem"] = str(config['resources']['gatk']['gatk_java_mem'])
    return

# old style bash script dispatch (for testing ==> ported into a ruffus based pipeline)
def molpathClassic(p):
    print >> sys.stderr, "DISPATCHING", p

    # create working directory
    syscall('mkdir',p.wd())

    # create SGE instance for mapping ()
    sge = SGE( out=None, err=None, queue=None, email=None, sendemail='a')

    # ALIGNMENT
    sge.setResources(multithread=4, vmem=8)

    # pipeline (auto-adds job dependencies / work only for linear pipeline)
    mp.trimAdapters(p,sge)
    mp.mapReads(p,sge)


    sge.setResources(multithread=4, vmem=8)
    mp.sam2bam(p,sge)
    mp.sortBam(p,sge)



    ## addreplaceReadgroups

    ## mark duplicates

    # GATK
    ## RealignerTargetCreator
    ## IndelRealigner
    ## BaseRecalibrator PRE-recalibration
    ## PrintReads (recalibration)
    ## BaseRecalibrator POST-recalibration

    # stampy/platypus

    #


    return



# main
if __name__=="__main__":
    # read requirements/includes
    usage = "usage: %prog [options] <sample descriptions>"
    parser = OptionParser(usage=usage)
    parser.add_option("-p", dest="pipepar", default='ngs_config.json',metavar="STRING", help="Pipeline parameters (ngs_config.json)")
    parser.add_option("-l", dest="logfile", default=None,             metavar="STRING", help="logfile (default STDERR)")
    (options, args) = parser.parse_args()
    if len(args) == 0:
        parser.print_help()
        sys.exit("\nERROR: no sample file provided")

    # start logging instance
    logging.basicConfig(
            filename=options.logfile,
            filemode='w',
            format='%(asctime)s - %(levelname)s - %(message)s',
            datefmt='%m/%d/%Y %I:%M:%S %p',
            level=logging.DEBUG)

    # read pipeline parameters
    logging.info("reading configuration from %s", options.pipepar)
    with open(options.pipepar) as pp:
        setupPipeline(json.load(pp))

    # process patient files
    for patientfile in args:
        # read sample/patient data
        with open(patientfile) as fh:
            patients = Patients(fh)

        # print patient data in readable format
        for i, patient in enumerate(patients):
            print i, patient
            print repr(patient)

            # setup pipeline jobs
            molpathClassic(patient)





