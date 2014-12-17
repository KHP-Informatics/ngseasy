


logger_ngseasy
logger_ngseasy_new
ngseasy

ngseasy_alignment_v0.9




ngseasy_functions



ngseasy_qcfiler_bam
ngseasy_filter_recalbam


ngs_full_gatk
ngs_full_no_gatk



## Step 1

## get contianers 
get_containers.sh

## step up directory structure 

ngseasy_initiate_project 

## copy fastq file from storage folder to project and sample folders
ngseasy_initiate_fastq

## start volumes container
ngseasy_volumes_container

### TO DO 
add checks for containers/images 

****

## NGS pipeline

```
ngseasy_fastqc

if [ trimm==1 ]
  then
    ngseasy_trimmomatic
fi

ngseasy_alignment

# GATK Processing
if [ gatk==1 ] && [ realign==1 ] 
  then 
    ngseasy_indel_realn
    ngseasy_base_recal
fi

## NON-GATK Processing
if [ gatk==0 ] && [ realign==1 ] 
  then 
    ngseasy_ogap_realn
    ngseasy_bamutil_base_recal
    
  else
    ngseasy_bamutil_base_recal
fi



ngseasy_alignment_qc

ngseasy_variant_calling

ngseasy_variant_calling_fast_ensemble


```


****

## Dumped

```
ngseasy_addreadgroup
ngseasy_markduplicates
```









