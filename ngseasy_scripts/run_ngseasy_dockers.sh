#!/bin/bash
config_csv=${1}
sudo docker run -d -P \
-v ../fastq_raw:/home/pipeman/fastq_raw \
-v ../reference_geneomes:/usr/local/pipeman/reference_genomes \
-v ../gatk_resources:/usr/local/pipeman/gatk_resources \
-v ../ngs_projects:/home/pipeman/ngs_projects \
-v ../ngseay_scripts:/usr/local/pipeman/ngseasy_scripts \
-u pipeman \
-t snewhouse/alignment-public:v1.2 /sbin/my_init -- bash /usr/local/pipeman/ngseasy_scripts/run-ea-ngs.sh /home/pipeman/ngs_projects/${config_csv}
