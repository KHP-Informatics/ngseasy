#########################################################################
# -- Author: Amos Folarin                                               #
# -- Email: amosfolarin@gmail.com                                       #
# -- Author: Stephen J Newhouse                                         #
# -- Email: stephen.j.newhouse@gmail.com                                #
# -- Organisation: KCL/SLaM                                             #
#########################################################################


#------------------------------------------------------------------------
# A simple script to acquire all the container images from DockerHub for
#  a given pipeline version.
#
# Either run as sudo or if you have docker configured to run without sudo
# then run as the regular user
# 
# Each major pipeline release has been given a tag, this file will enable 
# a record of which containers constituted the pipeline at that version
#------------------------------------------------------------------------

cat << EOF
# A simple script to acquire NGSeasy container images from DockerHub.
# Each major pipeline release has been given a tag, this file will enable
# a record of which containers constituted the pipeline at that version
# (see version tag in GitHub  https://github.com/KHP-Informatics/ngs )
#
# USAGE:
# Either run as sudo or if you have docker configured to run without sudo
# then run as the regular user and it will download the correct set of 
# container images
EOF

docker pull compbio/ngseasy-alignment-public:v1.2
docker pull compbio/ngseasy-platypus:v1.0
docker pull compbio/ngseasy-freebayes:v1.0
docker pull compbio/ngseasy-vcftools:v1.0
docker pull compbio/ngseasy-stampy:v1.0
docker pull compbio/ngseasy-samtools:v1.0
docker pull compbio/ngseasy-picard:v1.0
docker pull compbio/ngseasy-gatk:v1.0
docker pull compbio/ngseasy-bwa:v1.0
docker pull compbio/ngseasy-bowtie2:v1.0
docker pull compbio/ngseasy-bedtools:v1.0
docker pull compbio/ngseasy-sv-exomdepth:v1.0
docker pull compbio/ngseasy-sv-cnmops:v1.0
docker pull compbio/ngseasy-sv-mhmm:v1.0
docker pull compbio/ngseasy-sv-delly:v1.0
docker pull compbio/ngseasy-var-anno:v1.1
docker pull compbio/ngseasy-alignment-public:v1.1
docker pull compbio/ngseasy-fastqc:v1.0
docker pull compbio/ngseasy-base:v1.0


