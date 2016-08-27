#!/bin/bash

## Taken from Heng Li awesomeness
## https://github.com/lh3/bwa/blob/master/bwakit/run-gen-ref
root=`dirname $0`

# 1000genomes b38
url381kG="ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa" # 16.07.2016
url38="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.fna.gz"

# 1000genomes b37
url37d5="ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz"
url37="ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz"

# hg19 ucsc genome version for some tools
urlHG19="http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/chromFa.tar.gz"

if [ $# -eq 0 ]; then
	echo "Usage: $0 <hs38|hs38a|hs38DH|hs37|hs37d5|hg19>"
	echo "Analysis sets:"
	echo "  hs38     primary assembly of GRCh38 (incl. chromosomes, unplaced and unlocalized contigs) and EBV"
	echo "  hs38a    hs38 plus ALT contigs"
	echo "  hs38DH   hs38a plus decoy contigs and HLA genes (recommended for GRCh38 mapping)"
	echo "  hs37     primary assembly of GRCh37 (used by 1000g phase 1) plus the EBV genome"
	echo "  hs37d5   hs37 plus decoy contigs (used by 1000g phase 3)"
	echo ""
	echo "Note: This script downloads human reference genomes. For hs38a and hs38DH, it needs additional"
	echo "      sequences and ALT-to-REF mapping included in the bwa.kit package."
	exit 1;
fi

if [ $1 == "hs38DH" ]; then
	wget -O- $url38 | gzip -dc; cat $root/resource-GRCh38/hs38DH-extra.fa > $1.fa
	[ ! -f $1.fa.alt ] && cp $root/resource-GRCh38/hs38DH.fa.alt $1.fa.alt
elif [ $1 == "hs38a" ]; then
	wget -O- $url38 | gzip -dc > $1.fa
	[ ! -f $1.fa.alt ] && grep _alt $root/resource-GRCh38/hs38DH.fa.alt > $1.fa.alt
elif [ $1 == "hs38" ]; then
	wget -O- $url38 | gzip -dc | awk '/^>/{f=/_alt/?0:1}f' > $1.fa
elif [ $1 == "hs37d5" ]; then
	wget -O- $url37d5 | gzip -dc > $1.fa 2>/dev/null
elif [ $1 == "hs37" ]; then
	wget -O- $url37d5 | gzip -dc 2>/dev/null | awk '/^>/{f=/>hs37d5/?0:1}f' > $1.fa
elif [ $1 == "hg19" ]; then
  wget $urlHG19 -O $root/hg19/${1}.tar.gz
  tar xvzf $root/hg19/${1}.tar.gz
  for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M; do
    cat $root/hg19/chr${i}.fa >> $root/hg19/${1}.fa
    mv $root/hg19/chr${i}.fa $root/hg19/chromosomes/
  done
else
	echo "ERROR: unknown genome build"
fi


chr4_ctg9_hap1
chr6_apd_hap1
chr6_cox_hap2
chr6_dbb_hap3
chr6_mann_hap4
chr6_mcf_hap5
chr6_qbl_hap6
chr6_ssto_hap7
chr17_ctg5_hap1
chr1_gl000191_random
chr1_gl000192_random
chr4_gl000193_random
chr4_gl000194_random
chr7_gl000195_random
chr8_gl000196_random
chr8_gl000197_random
chr9_gl000198_random
chr9_gl000199_random
chr9_gl000200_random
chr9_gl000201_random
chr11_gl000202_random
chr17_gl000203_random
chr17_gl000204_random
chr17_gl000205_random
chr17_gl000206_random
chr18_gl000207_random
chr19_gl000208_random
chr19_gl000209_random
chr21_gl000210_random
chrUn_gl000211
chrUn_gl000212
chrUn_gl000213
chrUn_gl000214
chrUn_gl000215
chrUn_gl000216
chrUn_gl000217
chrUn_gl000218
chrUn_gl000219
chrUn_gl000220
chrUn_gl000221
chrUn_gl000222
chrUn_gl000223
chrUn_gl000224
chrUn_gl000225
chrUn_gl000226
chrUn_gl000227
chrUn_gl000228
chrUn_gl000229
chrUn_gl000230
chrUn_gl000231
chrUn_gl000232
chrUn_gl000233
chrUn_gl000234
chrUn_gl000235
chrUn_gl000236
chrUn_gl000237
chrUn_gl000238
chrUn_gl000239
chrUn_gl000240
chrUn_gl000241
chrUn_gl000242
chrUn_gl000243
chrUn_gl000244
chrUn_gl000245
chrUn_gl000246
chrUn_gl000247
chrUn_gl000248
chrUn_gl000249
