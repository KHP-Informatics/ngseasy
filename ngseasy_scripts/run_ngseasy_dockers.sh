#!/bin/bash

config_tsv=${1} 

sudo docker run -P \
-v /media/D/docker_ngs/ngseasy/fastq_raw:/home/pipeman/fastq_raw \
-v /media/D/docker_ngs/ngseasy/reference_genomes_b37:/home/pipeman/reference_genomes_b37 \
-v /media/D/docker_ngs/ngseasy/gatk_resources:/home/pipeman/gatk_resources \
-v /media/D/docker_ngs/ngseasy/ngs_projects:/home/pipeman/ngs_projects \
-v /media/D/docker_ngs/ngseasy/ngseasy_scripts:/home/pipeman/ngseasy_scripts \
-i -t snewhouse/ngseasy-alignment-public:v1.2 /sbin/my_init -- /bin/bash  /home/pipeman/ngseasy_scripts/run_ea-ngs.sh /home/pipeman/ngs_projects/${config_tsv}; chmod -R 755 /home/pipeman/;


