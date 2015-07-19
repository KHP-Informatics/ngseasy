#!/bin/bash


#################################################################################
# usage
#################################################################################

function usage {
    echo "
Program: VcfClean
Version 1.0-r001
Author: Stephen Newhouse (stephen.j.newhouse@gmail.com)

usage:   VcfClean -v <VCF> -f <Variant Caller>

options:  -v  STRING  configuration file
          -f  STRING  project directory
          -h  NULL    show this message
"
}

#################################################################################
# get options
#################################################################################

  while  getopts "hv:f:" opt
  do

      case ${opt} in
    h)
    usage #print help
    exit 0
    ;;

    c)
    VCF_INPUT=${OPTARG}
    ;;

    d)
    VCF_FORMAT=${OPTARG}
    ;;
      esac
  done


VCF_FORMAT="platypus"
VCF_INPUT="/media/Data/ngs_projects/GCAT_Data/NA12878/vcf/platypus_180715/small.vcf"

#################################################################################
## make compatible with other tools and mathc GATK/Freebayes annotaions
#################################################################################

if [[ ${VCF_FORMAT} == "platypus" ]]
    then
    if [[ gzip ]]
        then
        VCF=`basname ${VCF} .gz`
        zcat ${VCF}.gz | sed s/TC\=/DP\=/g | \
        sed s/"GT:GL:GOF:GQ:NR:NV"/"GT:GL:GOF:GQ:DP:NV"/g | bgzip -c > ${VCF}.fix.vcf.gz && \
        tabix ${VCF}.fix.vcf.gz
    elif [[ not gzip ]]
        then
        VCF=`basname ${VCF} .vcf`
        cat ${VCF}.vcf | sed s/TC\=/DP\=/g | \
        sed s/"GT:GL:GOF:GQ:NR:NV"/"GT:GL:GOF:GQ:DP:NV"/g |  bgzip -c > ${VCF}.fix.vcf.gz && \
        tabix ${VCF}.fix.vcf.gz
    else
        echo "ok"
    fi
else
    echo "ok"
fi

## TO FIX:- Get GATK/Freebayes equivalents and check GEMINI requirements
##FORMAT=<ID=NV,Number=.,Type=Integer,Description="Number of reads containing variant in this sample">
##FORMAT=<ID=NR,Number=.,Type=Integer,Description="Number of reads covering variant location in this sample">
##INFO=<ID=TCR,Number=1,Type=Integer,Description="Total reverse strand coverage at this locus">
##INFO=<ID=TCF,Number=1,Type=Integer,Description="Total forward strand coverage at this locus">
##INFO=<ID=TC,Number=1,Type=Integer,Description="Total coverage at this locus">
##INFO=<ID=NR,Number=.,Type=Integer,Description="Total number of reverse reads containing this variant">
##INFO=<ID=FR,Number=.,Type=Float,Description="Estimated population frequency of variant">
