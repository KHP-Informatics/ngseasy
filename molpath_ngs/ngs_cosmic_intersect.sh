#########################################################################
# -- Author: Amos Folarin                                               #
# -- Organisation: KCL/SLaM                                             #
# -- Email: amosfolarin@gmail.com                                       #
#########################################################################


#------------------------------------------------------------------------
# DESC:
#This gets the intersection of variants which are found in the cancer 
# variation database Catalogue of Somatic Mutations in Cancer (COSMIC) 
# ftp://ftp.sanger.ac.uk/pub/CGP/cosmic/
# vcf files for the mutations 
# ftp://ngs.sanger.ac.uk/production/cosmic/CosmicCodingMuts_v66_20130725.vcf.gz
# 
#This will provide information on drug sensitivities of each mutation,
# which while not definitive (they are based on in-vitro screens) will be 
# useful to a clinician looking at the output of this pipeline, there is 
# also a tonne of other information on mutation hotspots etc.
#
#There are some slides here on the project 
# Cancer Genome Project ftp://ftp.sanger.ac.uk/pub/seminars/sequencing_beyond/Cancer_Genome_Project.pdf 
#------------------------------------------------------------------------

# USAGE:

# ARGS:
vcf_patient=${1}  #patient vcf file from genotype caller
vcf_cosmic_c=${2}   #cosmic download vcf coding variants
vcf_cosmic_nc=${3} #cosmic download vcf non-coding variants

#intersect 
bedtools intersect -wa -s -a ${vcf_patient} -b ${vcf_cosmic_c} > cosmic_i-c_${vcf_patient}

bedtools intersect -wa -s -a ${vcf_patient} -b ${vcf_cosmic_nc} > cosmic_i-nc_${vcf_patient}


