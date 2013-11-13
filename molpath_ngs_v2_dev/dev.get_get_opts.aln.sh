

## options

while getopts "f:o:k:q:c:w:r:" OPTION; do
	case "${OPTION}" in

		f)	echo "-f fastq prefix is ${OPTARG}"				fastq=${OPTARG};;
		o)	echo "-o sam output is ${OPTARG}"				out_sam=${OPTARG};;
		k)	echo "-k novoalign k stats out is ${OPTARG}"	novo_k_stat=${OPTARG};;
		q)	echo "-q quality scores are ${OPTARG}"			novo_qual=${OPTARG};;
		c)	echo "-c number of cpus"						novo_core=${OPTARG};;
		w)	echo "-w working directory is ${OPTARG}"		working_dir=${OPTARG};;
		r)	echo "-r reference genome is ${OPTARG}"			reference_genome=${OPTARG};;
	esac

done

qsub bash.sh \
-f ${fastq_dir}/${fastq_prefix} \
-o ${sample_dir}/${sample_name}.aln.sam \
-k ${sample_name}.novoalign.K.stats \
-q STDFQ \
-c 8 \
-w ${sample_dir} \
-r ${reference_genome_novoindex};

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

