#!/bin/bash
set -x

#variables 
#config_tsv="" 
#host_vol_dir=""

#usage printing func
usage()
{
cat << EOF
This script sets up the docker volumes and starts the sequence alignment container:
See NGSEasy containerized instructions.

ARGUMENTS:
-h      Flag: Show this help message
-c      Config pipeline file
-d      Base directory for (fastq_raw, reference_genomes_b37, gatk_resources, ngs_projects, ngseasy_scripts)
EXAMPLE USAGE:
run_ngseasy_dockers.sh -c myconfigfile.tsv -d /media/D/docker_ngs/ngseasy/


EOF
}


#get options for command line args
while  getopts "h:c:d:" opt
do

    case ${opt} in
        h)
        usage #print help
        exit 0
        ;;
        
        c)
        config_tsv=${OPTARG}
        ;;
        
        d)
        host_vol_dir=${OPTARG}
        ;;

    esac
done


echo "${host_vol_dir}"

echo "${config_tsv}"


#check exists.
if [[ ! -e ${host_vol_dir}/fastq_raw ]] 
then
	echo " ${host_vol_dir}/fastq_raw does not exist "
	usage;
	exit 1;
fi

if [[ ! -e ${host_vol_dir}/reference_genomes_b37 ]]
then
	echo "${host_vol_dir}/reference_genomes_b37 does not exist"
	usage
	exit 1
fi

if [[ ! -e ${host_vol_dir}/gatk_resources ]]
then
	echo "${host_vol_dir}/gatk_resources does not exist" 
	usage
	exit 1
fi

if [[ ! -e ${host_vol_dir}/ngs_projects ]] 
then
	echo "${host_vol_dir}/ngs_projects does not exist"
	usage
	exit 1
fi

if [[ ! -e ${host_vol_dir}/ngseasy_scripts ]]
then
	echo "${host_vol_dir}/ngseasy_scripts does not exist"
	usage
	exit 1
fi

echo " Run Docker Image "

# run docker 
sudo docker run -P \
         -v ${host_vol_dir}/ngseasy/fastq_raw:/home/pipeman/fastq_raw \
         -v ${host_vol_dir}/ngseasy/reference_genomes_b37:/home/pipeman/reference_genomes_b37 \
         -v ${host_vol_dir}/ngseasy/gatk_resources:/home/pipeman/gatk_resources \
         -v ${host_vol_dir}/ngseasy/ngs_projects:/home/pipeman/ngs_projects \
         -v ${host_vol_dir}/ngseasy/ngseasy_scripts:/home/pipeman/ngseasy_scripts \
         -i -t snewhouse/ngseasy-alignment-public:v1.2 /sbin/my_init -- /bin/bash /home/pipeman/ngseasy_scripts/run_ea-ngs.sh /home/pipeman/ngs_projects/${config_tsv}; 




