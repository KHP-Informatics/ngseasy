#!/bin/bash


#################################################################################
# get options
#################################################################################

VCF="/media/Data/ngs_projects/GCAT_Data/NA12878/vcf/platypus_180715/small.vcf"

#################################################################################
## make compatible with other tools and mathc GATK/Freebayes annotaions
#################################################################################

if [[ gzip ]]
    then
        zcat ${VCF} | sed s/TC\=/DP\=/g | \
        sed s/"GT:GL:GOF:GQ:NR:NV"/"GT:GL:GOF:GQ:DP:NV"/g | bgzip -c > ${VCF}.fix.gz && \
        tabix ${VCF}.fix.gz
elif [[ not gzip ]]
    then
        cat ${VCF} | sed s/TC\=/DP\=/g | \
        sed s/"GT:GL:GOF:GQ:NR:NV"/"GT:GL:GOF:GQ:DP:NV"/g |  bgzip -c > ${VCF}.fix.gz && \
        tabix ${VCF}.fix.gz
else
    ok
fi


## TO FIX:-
##FORMAT=<ID=NV,Number=.,Type=Integer,Description="Number of reads containing variant in this sample">
##FORMAT=<ID=NR,Number=.,Type=Integer,Description="Number of reads covering variant location in this sample">
##INFO=<ID=TCR,Number=1,Type=Integer,Description="Total reverse strand coverage at this locus">
##INFO=<ID=TCF,Number=1,Type=Integer,Description="Total forward strand coverage at this locus">
##INFO=<ID=TC,Number=1,Type=Integer,Description="Total coverage at this locus">
##INFO=<ID=NR,Number=.,Type=Integer,Description="Total number of reverse reads containing this variant">
##INFO=<ID=FR,Number=.,Type=Float,Description="Estimated population frequency of variant">
