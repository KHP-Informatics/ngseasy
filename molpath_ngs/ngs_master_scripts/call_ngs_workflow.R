#!/share/bin/Rscript --vanilla --default-packages=utils

args <- commandArgs(trailingOnly=TRUE);

config_file <- args[1];

d <- read.table(config_file, head=T,sep="\t",as.is=T,fill=T)


for( i in 1:dim(d)[1] )

{

pipeline <- d$Pipeline[i] ## pipeline to call
Fastq_prefix=d$Fastq_file_prefix[i]	## fastq prefix
sample_name=d$ReadGroup_sample_RGSM[i] ##sample name
qual_type=d$QUAL[i]  ## Base quality coding for novoalign ie STFQ, ILMFQ, ILM1.8
RGID=d$ReadGroup_id_RGID[i] 	#Read Group ID Required.
RGLB=d$ReadGroup_library_RGLB[i] 	#Read Group Library Required.
RGPL=d$ReadGroup_platform_RGPL[i] 	#Read Group platform (e.g. illumina, solid) Required.
RGPU=d$ReadGroup_platform_unit_RGPU[i] 	#Read Group platform unit (eg. run barcode) Required.
RGSM=d$ReadGroup_sample_RGSM[i] 	#Read Group sample name Required.
RGCN=d$ReadGroup_SeqCentre_RGCN[i] 	#Read Group sequencing center name Required.
RGDS=d$ReadGroup_Desc_RGDS[i] 	#Read Group description Required.
RGDT=d$ReadGroup_runDate_RGDT[i] 	#Read Group run date Required.
isPE=d$PE[i] 	#Read Group run date Required.
targetbed=d$bed_list[i]
bedtype=d$bed_type[i]
email_user=d$email_bioinf[i]
fastq_dir=d$Fastq_dir[i]
bamdir=d$BAM_dir[i]

callPipe <- paste(" qsub -N call_ngs.",RGID," ", pipeline," ",Fastq_prefix," ", sample_name," ",qual_type," ",RGID," ",RGLB," ",RGPL," ",RGPU," ",RGSM," ",RGCN," ",RGDS," ",RGDT," ",isPE," ",targetbed," ",bedtype," ",email_user," ",fastq_dir," ",bamdir,sep="")

system(callPipe)

cat(callPipe,"\r","\n")

}



