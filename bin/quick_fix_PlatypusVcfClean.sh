#!/bin/bash


#################################################################################
# get options
#################################################################################

VCF=""

#################################################################################
## make compatible with other tools and mathc GATK/Freebayes annotaions
#################################################################################

if [[ gzip ]]
    then
        zcat ${VCF} | sed s/TC/DP/g | sed s/"GT:GL:GOF:GQ:NR:NV"/"GT:GL:GOF:GQ:DP:NV"/g | bgzip -c > ${VCF}.fix.gz && \
        tabix ${VCF}.fix.gz
elif [[ not gzip ]]
    then
        cat ${VCF} | sed s/TC/DP/g | sed s/"GT:GL:GOF:GQ:NR:NV"/"GT:GL:GOF:GQ:DP:NV"/g |  bgzip -c > ${VCF}.fix.gz && \
        tabix ${VCF}.fix.gz
else
    ok
fi
