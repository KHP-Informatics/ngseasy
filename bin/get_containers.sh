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

VERSION=${1}

sudo docker pull compbio/ngseasy-base:latest
sudo docker pull compbio/ngseasy-fastqc:${VERSION}
sudo docker pull compbio/ngseasy-trimmomatic:${VERSION}
sudo docker pull compbio/ngseasy-bwa:${VERSION}
sudo docker pull compbio/ngseasy-bowtie:${VERSION}
sudo docker pull compbio/ngseasy-picardtools:${VERSION}
sudo docker pull compbio/ngseasy-samtools:${VERSION}
sudo docker pull compbio/ngseasy-freebayes:${VERSION}
sudo docker pull compbio/ngseasy-bedtools:${VERSION}
sudo docker pull compbio/ngseasy-bcbiovar:${VERSION}
sudo docker pull compbio/ngseasy-delly:${VERSION}
sudo docker pull compbio/ngseasy-lumpy:${VERSION}
sudo docker pull compbio/ngseasy-cnmops:${VERSION}
sudo docker pull compbio/ngseasy-mhmm:${VERSION}
sudo docker pull compbio/ngseasy-exomedepth:${VERSION}
sudo docker pull compbio/ngseasy-bamutil:${VERSION}


