#!/bin/bash -x

RUNDATE="230715"

NGSEASY_DIR="/media/Data/ngs_projects"

for aln in bwa snap stampy novoalign bowtie2; do
    for qctrim in no-trim; do
        for vcaller in platypus freebayes-parallel; do

echo -e "\n\n[`date`]Start NGSeasy : $aln $qctim $vcaller \n\n "
	
	NGSEASY_CONFIG="${aln}.${qctrim}.${vcaller}.${RUNDATE}.config.tsv"
	
	LOGFILE="${aln}.${qctrim}.${vcaller}.${RUNDATE}.run.log"
        touch ${LOGFILE}
        
	NGSEASY_EXEC=`which ngseasy`

	time ${NGSEASY_EXEC} -c ${NGSEASY_CONFIG} -d ${NGSEASY_DIR} | tee -a >> ${LOGFILE};
        
	sleep 1s
        
	rm ./*tmp

echo -e "\n\n[`date`]END NGSeasy : $aln $qctim $vcaller \n\n "
        
	done
    done
done

