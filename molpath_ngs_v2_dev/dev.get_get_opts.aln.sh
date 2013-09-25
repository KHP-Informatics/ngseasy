

## options

while getopts fq:os:ok:qu:nc:wd:rg: option; do
	case "${option}" in

		fq)		fastq=${OPTARG};;
		os)		out_sam=${OPTARG};;
		ok)		novo_k_stat=${OPTARG};;
		qu)		novo_qual=${OPTARG};;
		nc)		novo_core=${OPTARG};;
		nc)		working_dir=${OPTARG};;
		rg)		reference_genome=${OPTARG};;
	esac

done

qsub bash.sh \
-fq ${fastq_dir}/${fastq_prefix} \
-os ${sample_dir}/${sample_name}.aln.sam \
-ok ${sample_name}.novoalign.K.stats \
-qu STDFQ \
-nc 8 \
-wd ${sample_dir} \
-wd ${reference_genome_novoindex};

#---------------------------------------------------------------------------------------------------#

cd ${working_dir}

${ngs_novo}/novoalign \
-d ${reference_genome} \
-f ${fastq}_1.fastq  ${fastq}_2.fastq  \
-F ${novo_qual} \
--Q2Off \
--3Prime  \
-g 40 \
-x 6 \
-r All \
-i PE 300,150 \
-a ACACTCTTTCCCTACACGACGCTCTTCCGATCT GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG \
-a AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT  CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT \
-a ACACTCTTTCCCTACACGACGCTCTTCCGATCT CGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT \
-a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA CGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT \
-a AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA \
-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
-c ${novo_core} \
-k -K ${novo_k_stat} \
-o SAM > ${out_sam};

echo "end novoalign"

echo "----------------------------------------------------------------------------------------"

