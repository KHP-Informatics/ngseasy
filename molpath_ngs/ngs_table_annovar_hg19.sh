#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -p -0.99999999999999999999999999999999999999999999999999
#$ -V

sample_name=${1}
sample_dir=${2}
sample_temp=${3}
geno=${4}

cd ${sample_dir}

invcf=${sample_dir}/${sample_name}.novorecal.${geno}.raw.snps.indels.vcf
avinput=${sample_dir}/${sample_name}.novorecal.${geno}.raw.snps.indels.avinput
table_annovar_out=${sample_dir}/${sample_name}.novorecal.${geno}.raw.snps.indels

##################################
## convert vcf to annovar input ##
##################################

${annovar}/convert2annovar.pl --includeinfo --genoqual 0 --varqual 0 --format vcf4 ${invcf} --outfile ${avinput};

##                --withzyg                   print zygosity/coverage/quality when -includeinfo is used (for vcf4 format)
##                --genoqual <float>          genotype quality score threshold (for vcf4 format)
##                --varqual <float>           variant quality score threshold (for vcf4 format)
##                --comment                   keep comment line in output (for vcf4 format)



#######################
## run table annovar ##
#######################

${annovar}/table_annovar.pl ${avinput} ${annovar_humandb} \
--protocol refGene,knownGene,ensGene,wgEncodeGencodeManualV4,gerp++elem,phastConsElements46way,genomicSuperDups,esp6500si_all,1000g2012apr_all,1000g2012apr_eur,1000g2012apr_amr,1000g2012apr_asn,1000g2012apr_afr,cg46,cosmic64,snp129,snp132,snp137,avsift,ljb2_all \
--operation g,g,g,g,r,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f \
--csvout \
--otherinfo \
--buildver hg19 \
--remove \
--outfile ${table_annovar_out}; 




#Since June 2013, the LJB2_* databases are now provided to ANNOVAR users. 
## Therefore, users can try replace the ljb_all in the command above with ljb2_all, and will get more detailed exome annotations. 
## In addition, the -arg argument is now supported, 
## so that you can supply a list of comma-delimited optional arguments to table_annovar 
## for each of the annotation tasks. 
## For example, adding -arg '-splicing 5',,,,,,,,,,,, to the command will add change the splicing threshold to 5bp for the gene-based annotation.


