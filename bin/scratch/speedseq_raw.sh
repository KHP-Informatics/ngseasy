#!/bin/bash -e

############################################################
#  Program: speedseq
#  Version: 0.0.2
#  Author: Colby Chiang (cc2qe@virginia.edu)
############################################################

# source the paths to the binaries used in the script
function source_binaries() {
    if [[ -e $1 ]]
    then
	echo "Sourcing executables from $1 ..."
	if [[ $1 == /* ]]
	then
	    source $1
	else
	    source ./$1
	fi
    else
	echo "Config file $1 not found. Attempting to auto-source executables"
	# general
	SPEEDSEQ_HOME=$( dirname `which speedseq` )
	SAMBAMBA=`which sambamba || true`
	BEDTOOLS=`which bedtools || true`
	BGZIP=`which bgzip || true`
	TABIX=`which tabix || true`
	VAWK=`which vawk || true`
	PARALLEL=`which parallel || true`
	PYTHON27=`which python2.7 || true`

	# align
	BWA=`which bwa || true`
	SAMBLASTER=`which samblaster || true`

	# var/somatic
	FREEBAYES=`which freebayes || true`
	VEP=`which variant_effect_predictor.pl || true`
	VEP_CACHE_DIR=$SPEEDSEQ_HOME/annotations/vep_cache

	# sv
	LUMPY=`which lumpy || true`
	PAIREND_DISTRO=`which pairend_distro.py || true`
	BEDPETOVCF=`which bedpeToVcf || true`
	LUMPYTOBEDPE=`which lumpyToBedpe || true`
	SVTYPER=`which svtyper || true`

        # CNVnator
	CNVNATOR_WRAPPER=`which cnvnator_wrapper.py || true`
	CNVNATOR_MULTI=`which cnvnator-multi || true`
	ANNOTATE_RD=`which annotate_rd.py || true`
	CNVNATOR_CHROMS_DIR=~/genomes/GRCh37/chroms

	# re-align
	BAMTOFASTQ=`which bamtofastq.py || true`
	MBUFFER=`which mbuffer || true`
	BAMHEADRG=`which bamheadrg.py || true`
    fi
}

# ensure that the require python modules are installed before
# beginning analysis
function check_python_modules() {
    PYTHON_TEST=$1
    echo -e "\nChecking for required python modules ($PYTHON_TEST)..."

    $PYTHON_TEST -c "import imp; imp.find_module('pysam')"
    $PYTHON_TEST -c "import imp; imp.find_module('numpy')"
    $PYTHON_TEST -c "import imp; imp.find_module('scipy')"
}

## global usage
function usage() {
    echo "
Program: speedseq
Version: 0.0.1
Author: Colby Chiang (cc2qe@virginia.edu)

usage:   speedseq <command> [options]

command: align    align FASTQ files with BWA-MEM
         var      call SNV and indel variants with FreeBayes
         somatic  call somatic SNV and indel variants in a tumor/normal pair with FreeBayes
         sv       call SVs with LUMPY
         realign  re-align from a coordinate sorted BAM file

options: -h       show this message
"
}

function somatic_filter() {
    awk -v MINQUAL="$1" -v SSC_THRES="$2" -v ONLY_SOMATIC="$3" 'BEGIN {NORMAL=10; TUMOR=11; GL_IDX=0;}
    {
        if ($0~"^#") { print ; next; }
        if (! GL_IDX) {
            split($9,fmt,":")
            for (i=1;i<=length(fmt);++i) { if (fmt[i]=="GL") GL_IDX=i }
        }
        split($NORMAL,N,":");
        split(N[GL_IDX],NGL,",");
        split($TUMOR,T,":");
        split(T[GL_IDX],TGL,",");
        LOD_NORM=NGL[1]-NGL[2];
        LOD_TUMOR_HET=TGL[2]-TGL[1];
        LOD_TUMOR_HOM=TGL[3]-TGL[1];

        if (LOD_TUMOR_HET > LOD_TUMOR_HOM) { LOD_TUMOR=LOD_TUMOR_HET }
        else { LOD_TUMOR=LOD_TUMOR_HOM }

        DQUAL=LOD_TUMOR+LOD_NORM;

        if (DQUAL>=SSC_THRES && $NORMAL~"^0/0") {
            $7="PASS"
            $8="SSC="DQUAL";"$8
            print
        }
        else if (!ONLY_SOMATIC && $6>=MINQUAL && $10~"^0/0" && ! match($11,"^0/0")) {
            $8="SSC="DQUAL";"$8
            print
        }
    }' OFS="\t"
}

# alignment with BWA-MEM
function align() {
    function align_usage() {
	echo "
usage:   speedseq align [options] <reference.fa> <in1.fq> [in2.fq]

positional args:
         reference.fa
                  fasta file (indexed with bwa)
         in1.fq   paired-end fastq file. if -p flag is used then expected to be
                    an interleaved paired-end fastq file, and in2.fq may be omitted.
                    (can be gzipped)
         in2.fq   paired-end fastq file. (can be gzipped)

alignment options:
         -o STR   output prefix [in1.fq]
         -R       read group header line such as \"@RG\tID:libraryname\tSM:samplename\" (required)
         -p       first fastq file consists of interleaved paired-end sequences
         -t INT   threads [1]
         -T DIR   temp directory [./output_prefix.XXXXXXXXXXXX]
         -I FLOAT[,FLOAT[,INT[,INT]]]
                  specify the mean, standard deviation (10% of the mean if absent), max
                    (4 sigma from the mean if absent) and min of the insert size distribution.
                    FR orientation only. [inferred]

samblaster options:
         -i       include duplicates in splitters and discordants
         -c INT   maximum number of split alignments for a read to be included in splitter file [2]
         -m INT   minimum non-overlapping base pairs between two alignments for a read to be included in splitter file [20]

sambamba options:
         -M       amount of memory in GB to be used for sorting [20]

global options:
         -K FILE  path to speedseq.config file (default: same directory as speedseq)
         -v       verbose
         -h       show this message
"
    }

    # Check options passed in.
    if test -z "$2"
    then
	align_usage
	exit 1
    fi

    # set defaults
    SPEEDSEQ_DIR=`dirname $0`
    CONFIG="$SPEEDSEQ_DIR/speedseq.config"
    INTERLEAVED=0
    RG_FMT=""
    OUTPUT=""
    INCLUDE_DUPS="--excludeDups"
    MAX_SPLIT_COUNT=2
    MIN_NON_OVERLAP=20
    THREADS=1
    TEMP_DIR=""
    VERBOSE=0
    INS_DIST=""
    SORT_MEM=20 # amount of memory for sorting, in gigabytes

    while getopts ":hw:o:R:pic:m:M:t:T:I:vK:" OPTION
    do
	case "${OPTION}" in
	    h)
		align_usage
		exit 1
		;;
	    R)
		RG="$OPTARG"
		RG_FMT="-R '$OPTARG'"
		;;
	    p)
		INTERLEAVED=1
		;;
	    o)
		OUTPUT="$OPTARG"
		;;
	    i)
		INCLUDE_DUPS=""
		;;
	    c)
		MAX_SPLIT_COUNT="$OPTARG"
		;;
	    m)
		MIN_NON_OVERLAP="$OPTARG"
		;;
	    M)
		SORT_MEM="$OPTARG"
		;;
	    t)
		THREADS="$OPTARG"
		;;
	    T)
		TEMP_DIR="$OPTARG"
		;;
	    I)
		INS_DIST="-I $OPTARG"
		;;
	    v)
		VERBOSE=1
		;;
	    K)
		CONFIG="$OPTARG"
		;;
	esac
    done

    if [[ "$INTERLEAVED" -eq 1 ]]
    then
	REF="${@:${OPTIND}:1}"
	FQ="${@:$((${OPTIND}+1)):1}"
	if [[ -z "$OUTPUT" ]]
	then
	    OUTPUT=`basename "$FQ"`
	fi

        # Check that the ref and fastq files exist
	if [[ -z "$REF" ]] || [[ ! -f "$REF" ]]
	then
            align_usage
            echo -e "Error: Reference file $REF not found.\n"
            exit 1
	elif [[ -z "$FQ" ]] || [[ ! -e "$FQ" ]]
	then
            align_usage
            echo -e "Error: Fastq file $FQ not found.\n"
            exit 1
	fi
    else
	REF="${@:${OPTIND}:1}"
	FQ1="${@:$((${OPTIND}+1)):1}"
	FQ2="${@:$((${OPTIND}+2)):1}"
	if [[ -z "$OUTPUT" ]]
	then
	    OUTPUT=`basename "$FQ1"`
	fi

        # Check that the ref and fastq files exist
	if [[ -z "$REF" ]] || [[ ! -f "$REF" ]]
	then
            align_usage
            echo -e "Error: Reference file $REF not found.\n"
            exit 1
	elif [[ -z "$FQ1" ]] || [[ ! -e "$FQ1" ]]
	then
            align_usage
            echo -e "Error: Fastq file $FQ1 not found.\n"
            exit 1
	elif [[ -z "$FQ2" ]] || [[ ! -e "$FQ2" ]]
	then
            align_usage
            echo -e "Error: Fastq file $FQ2 not found. (single-end reads not supported)\n"
            exit 1
	fi
    fi

    OUTBASE=`basename "$OUTPUT"`

    # Check for readgroup flag
    if [[ -z $RG_FMT ]]
    then
	align_usage
	echo -e "Error: no readgroup found. Please set a readgroup with the -R flag.\n"
	exit 1
    fi

    # Check the for the relevant binaries
    source_binaries $CONFIG

    if [[ -z "$BWA" ]]
    then
	align_usage
        echo -e "Error: bwa executable not found. Please set path in $SPEEDSEQ_DIR/speedseq.config file\n"
        exit 1
    elif [[ -z  "$SAMBLASTER" ]]
    then
	align_usage
        echo -e "Error: samblaster executable not found. Please set path in $SPEEDSEQ_DIR/speedseq.config file\n"
        exit 1
    elif [[ -z "$SAMBAMBA" ]]
    then
	align_usage
        echo -e "Error: sambamba executable not found. Please set path in $SPEEDSEQ_DIR/speedseq.config file\n"
        exit 1
    elif [[ -z "$PARALLEL" ]]
    then
	align_usage
        echo -e "Error: parallel executable not found. Please set path in $SPEEDSEQ_DIR/speedseq.config file\n"
        exit 1
    fi

    echo "Aligning..."
    # create temp directory if not specified by command argument
    if [[ -z $TEMP_DIR ]]
    then
	TEMP_DIR=`mktemp -d -p ./ ${OUTBASE}.XXXXXXXXXXXX`
    fi

    if [[ $VERBOSE -eq 1 ]]
    then
	echo "
        mkdir -p $TEMP_DIR/full $TEMP_DIR/spl $TEMP_DIR/disc
        mkfifo $TEMP_DIR/spl_pipe $TEMP_DIR/disc_pipe"
    fi

    # create temp files
    mkdir -p $TEMP_DIR/full $TEMP_DIR/spl $TEMP_DIR/disc
    if [[ ! -e $TEMP_DIR/spl_pipe ]]
    then
	mkfifo $TEMP_DIR/spl_pipe
    fi
    if [[ ! -e $TEMP_DIR/disc_pipe ]]
    then
	mkfifo $TEMP_DIR/disc_pipe
    fi

    # alignment command
    if [[ "$INTERLEAVED" -eq 1 ]]
    then
	if [[ $VERBOSE -eq 1 ]]
	then
	    echo -e "
        $BWA mem -t $THREADS -M -p $INS_DIST $RG_FMT $REF $FQ | \\
            $SAMBLASTER $INCLUDE_DUPS --addMateTags --maxSplitCount $MAX_SPLIT_COUNT --minNonOverlap $MIN_NON_OVERLAP --splitterFile $TEMP_DIR/spl_pipe --discordantFile $TEMP_DIR/disc_pipe | \\
            $SAMBAMBA view -S -f bam -l 0 /dev/stdin | \\
            $SAMBAMBA sort -t $THREADS -m $((${SORT_MEM}-2))G --tmpdir=$TEMP_DIR/full -o $OUTPUT.bam /dev/stdin

        $SAMBAMBA view -S -f bam -l 0 $TEMP_DIR/spl_pipe | \\
            $SAMBAMBA sort -t 4 -m 1G --tmpdir=$TEMP_DIR/spl -o $OUTPUT.splitters.bam /dev/stdin
        $SAMBAMBA view -S -f bam $TEMP_DIR/disc_pipe | \\
            $SAMBAMBA sort -t 4 -m 1G --tmpdir=$TEMP_DIR/disc -o $OUTPUT.discordants.bam /dev/stdin"
	fi

	echo "
        $BWA mem -t $THREADS -M -p $INS_DIST $RG_FMT $REF $FQ | \
	    $SAMBLASTER $INCLUDE_DUPS --addMateTags --maxSplitCount $MAX_SPLIT_COUNT --minNonOverlap $MIN_NON_OVERLAP --splitterFile $TEMP_DIR/spl_pipe --discordantFile $TEMP_DIR/disc_pipe | \
            $SAMBAMBA view -S -f bam -l 0 /dev/stdin | \
	    $SAMBAMBA sort -t $THREADS -m $((${SORT_MEM}-2))G --tmpdir=$TEMP_DIR/full -o $OUTPUT.bam /dev/stdin
        
        $SAMBAMBA view -S -f bam -l 0 $TEMP_DIR/spl_pipe | \
	    $SAMBAMBA sort -t 4 -m 1G --tmpdir=$TEMP_DIR/spl -o $OUTPUT.splitters.bam /dev/stdin
        $SAMBAMBA view -S -f bam $TEMP_DIR/disc_pipe | \
	    $SAMBAMBA sort -t 4 -m 1G --tmpdir=$TEMP_DIR/disc -o $OUTPUT.discordants.bam /dev/stdin
        " | $PARALLEL -j 3
    else
	if [[ $VERBOSE -eq 1 ]]
	then
        echo -e "
        $BWA mem -t $THREADS -M $INS_DIST $RG_FMT $REF $FQ1 $FQ2 | \\
            $SAMBLASTER $INCLUDE_DUPS --addMateTags --maxSplitCount $MAX_SPLIT_COUNT --minNonOverlap $MIN_NON_OVERLAP --splitterFile $TEMP_DIR/spl_pipe --discordantFile $TEMP_DIR/disc_pipe | \\
            $SAMBAMBA view -S -f bam -l 0 /dev/stdin | \\
            $SAMBAMBA sort -t $THREADS -m $((${SORT_MEM}-2))G --tmpdir=$TEMP_DIR/full -o $OUTPUT.bam /dev/stdin

        $SAMBAMBA view -S -f bam -l 0 $TEMP_DIR/spl_pipe | \\
            $SAMBAMBA sort -t 4 -m 1G --tmpdir=$TEMP_DIR/spl -o $OUTPUT.splitters.bam /dev/stdin
        $SAMBAMBA view -S -f bam $TEMP_DIR/disc_pipe | \\
            $SAMBAMBA sort -t 4 -m 1G --tmpdir=$TEMP_DIR/disc -o $OUTPUT.discordants.bam /dev/stdin"
	fi

        echo "
        $BWA mem -t $THREADS -M $INS_DIST $RG_FMT $REF $FQ1 $FQ2 | \
            $SAMBLASTER $INCLUDE_DUPS --addMateTags --maxSplitCount $MAX_SPLIT_COUNT --minNonOverlap $MIN_NON_OVERLAP --splitterFile $TEMP_DIR/spl_pipe --discordantFile $TEMP_DIR/disc_pipe | \
            $SAMBAMBA view -S -f bam -l 0 /dev/stdin | \
            $SAMBAMBA sort -t $THREADS -m $((${SORT_MEM}-2))G --tmpdir=$TEMP_DIR/full -o $OUTPUT.bam /dev/stdin

        $SAMBAMBA view -S -f bam -l 0 $TEMP_DIR/spl_pipe | \
            $SAMBAMBA sort -t 4 -m 1G --tmpdir=$TEMP_DIR/spl -o $OUTPUT.splitters.bam /dev/stdin
        $SAMBAMBA view -S -f bam $TEMP_DIR/disc_pipe | \
            $SAMBAMBA sort -t 4 -m 1G --tmpdir=$TEMP_DIR/disc -o $OUTPUT.discordants.bam /dev/stdin
        " | $PARALLEL -j 3
    fi
    
    # index the files
    if [[ $VERBOSE -eq 1 ]]
    then
	echo -e "
        $SAMBAMBA index $OUTPUT.bam
            $SAMBAMBA index $OUTPUT.discordants.bam
            $SAMBAMBA index $OUTPUT.splitters.bam"
    fi

    echo "
    $SAMBAMBA index $OUTPUT.bam
    $SAMBAMBA index $OUTPUT.discordants.bam
    $SAMBAMBA index $OUTPUT.splitters.bam
    " | $PARALLEL -j 3

    # clean up
    rm -r $TEMP_DIR

    echo "Done"

    # exit cleanly
    exit 0
}

function var() {
    function var_usage() {
        echo "
usage:   speedseq var [options] <reference.fa> <input1.bam> [input2.bam [...]]

positional args:
         reference.fa
                  genome reference fasta file
         input.bam
                  BAM file(s) to call variants on. Must have readgroup information,
                    and the SM readgroup tags will be the VCF column header

options:
         -o STR   output prefix [input1.bam]
         -w FILE  BED file of windowed genomic intervals
         -q FLOAT minimum variant QUAL score to output [1]
         -t INT   threads [1]
         -T DIR   temp directory [./output_prefix.XXXXXXXXXXXX]
         -A BOOL  annotate the vcf with VEP (true or false) (default: true)
         -K FILE  path to speedseq.config file (default: same directory as speedseq)
         -v       verbose
         -k       keep temporary files
         -h       show this message
"
    }

    # Check options passed in.
    if test -z "$2"
    then
        var_usage
        exit 1
    fi

    # set defaults
    SPEEDSEQ_DIR=`dirname $0`
    CONFIG="$SPEEDSEQ_DIR/speedseq.config"
    THREADS=1
    MINQUAL=1
    TEMP_DIR=""
    ANNOTATE="true"
    MINQUAL=1
    VERBOSE=0
    KEEP=0

    while getopts ":ho:w:t:T:A:q:vkK:" OPTION
    do
        case "${OPTION}" in
            h)
                var_usage
                exit 1
                ;;
            o)
                OUTPUT="$OPTARG"
                ;;
	    w)
		WINDOWS="$OPTARG"
		;;
	    q)
		MINQUAL="$OPTARG"
		;;
            t)
                THREADS="$OPTARG"
                ;;
            T)
                TEMP_DIR="$OPTARG"
                ;;
	    A)
		ANNOTATE=`echo "$OPTARG" | tr [:upper:] [:lower:]`
		;;
            v)
                VERBOSE=1
                ;;
	    k)
		KEEP=1
		;;
	    K)
		CONFIG="$OPTARG"
		;;
        esac
    done

    # parse the positional arguments
    REF="${@:${OPTIND}:1}"
    BAM_STRING="${@:$((${OPTIND}+1))}"
    BAM_LIST=($BAM_STRING)
    if [[ -z $OUTPUT ]]
    then
	OUTPUT=`basename "${BAM_LIST[0]}"`
    fi
    OUTBASE=`basename "$OUTPUT"`

    OPTIND=0

    # Check the for the relevant binaries
    source_binaries $CONFIG

    if [[ ! -f "$FREEBAYES" ]]
    then
	var_usage
        echo -e "Error: freebayes executable not found. Please set path in $SPEEDSEQ_DIR/speedseq.config file\n"
        exit 1
    elif [[ ! -f "$BEDTOOLS" ]]
    then
	var_usage
        echo -e "Error: bedtools executable not found. Please set path in $SPEEDSEQ_DIR/speedseq.config file\n"
        exit 1
    elif [[ ! -f "$BGZIP" ]]
    then
	var_usage
        echo -e "Error: bgzip executable not found. Please set path in $SPEEDSEQ_DIR/speedseq.config file\n"
        exit 1
    elif [[ ! -f "$TABIX" ]]
    then
	var_usage
        echo -e "Error: tabix executable not found. Please set path in $SPEEDSEQ_DIR/speedseq.config file\n"
        exit 1
    elif [[ ! -f "$VEP" ]] && [[ "$ANNOTATE" == "true" ]]
    then
	var_usage
        echo -e "Error: VEP not found. Please set path in $SPEEDSEQ_DIR/speedseq.config file\n"
        exit 1
    elif [[ ! -d "$VEP_CACHE_DIR" ]] && [[ "$ANNOTATE" == "true" ]]
    then
	var_usage
        echo -e "Error: VEP cache directory not found. Please set path in $SPEEDSEQ_DIR/speedseq.config file\n"
        exit 1
    elif [[ -z "$PARALLEL" ]]
    then
	var_usage
        echo -e "Error: parallel executable not found. Please set path in $SPEEDSEQ_DIR/speedseq.config file\n"
        exit 1
    elif [[ -z "$VAWK" ]]
    then
	var_usage
        echo -e "Error: vawk executable not found. Please set path in $SPEEDSEQ_DIR/speedseq.config file\n"
        exit 1
    fi

    # Check that ANNOTATE is either true or false
    if [[ "$ANNOTATE" != "true" ]] && [[ "$ANNOTATE" != "false" ]]
    then
	var_usage
	echo -e "Error: -A must be either true or false\n"
	exit 1
    fi

    # Check that the ref and bam files exist
    if [[ -z "$REF" ]] || [[ ! -f "$REF" ]]
    then
	var_usage
	echo -e "Error: Reference file $REF not found.\n"
	exit 1
    fi

    for TEST_BAM in ${BAM_LIST[@]}
    do
	if [[ ! -f $TEST_BAM ]]
	then
	    var_usage
	    echo -e "Error: BAM file $TEST_BAM not found.\n"
	    exit 1
	fi
    done

    echo "Calling variants..."
    # make temporary directory
    if [[ $VERBOSE -eq 1 ]]
    then
	echo "
        create temporary directory"
    fi
    if [[ -z $TEMP_DIR ]]
    then
	TEMP_DIR=`mktemp -d -p ./ ${OUTBASE}.XXXXXXXXXXXX`
    else
	mkdir -p $TEMP_DIR	
    fi

    # if no windows file, then make naive windows based on the chroms
    if [[ -z $WINDOWS ]]
    then
	if [[ $VERBOSE -eq 1 ]]
	then
	    echo -e "
        $SAMBAMBA view -H ${BAM_LIST[0]} | grep \"^@SQ\" | cut -f 2- | awk '{ gsub(\"^SN:\",\"\",\$1); gsub(\"^LN:\",\"\",\$2); print \$1\"\\\t0\\\t\"\$2; }' > $TEMP_DIR/windows.bed"
	fi

	$SAMBAMBA view -H ${BAM_LIST[0]} | grep "^@SQ" | cut -f 2- | awk '{ gsub("^SN:","",$1); gsub("^LN:","",$2); print $1"\t0\t"$2; }' > $TEMP_DIR/windows.bed
	WINDOWS="$TEMP_DIR/windows.bed"
    fi
    
    # construct the parallel command over windows
    if [[ $VERBOSE ]]
    then
        echo -e "
        $FREEBAYES \\
            -f $REF \\
            --region \$chrom:\$start..\$end \\
            $BAM_STRING \\
            --experimental-gls \\
            --min-repeat-entropy 1 \\
            | $VAWK --header '\$6>=$MINQUAL && I\$RPR>0 && I\$RPL>0' \\
            > ${TEMP_DIR}/$OUTBASE.\$i.vcf"
    fi
    
    for i in `cat $WINDOWS | awk '{print $1":"$2".."$3}'`
    do
        echo "$FREEBAYES \
        -f $REF \
        --region $i \
        --experimental-gls \
        --min-repeat-entropy 1 \
        $BAM_STRING \
        | $VAWK --header '\$6>=$MINQUAL && I\$RPR>0 && I\$RPL>0' \
        > ${TEMP_DIR}/$OUTBASE.$i.vcf"
    done > $TEMP_DIR/var_command.txt

    # run the parallel freebayes command
    if [[ $VERBOSE ]]
    then
	echo "
        cat $TEMP_DIR/var_command.txt | $PARALLEL -j $THREADS"
    fi
    cat $TEMP_DIR/var_command.txt | $PARALLEL -j $THREADS

    # make vcf header
    i=`head -n 1 $WINDOWS | awk '{print $1":"$2".."$3}'`
    if [[ $VERBOSE -eq 1 ]]
    then
	echo -e "
        grep \"^#\" $TEMP_DIR/$OUTBASE.$i.vcf > $TEMP_DIR/header.txt"
    fi
    grep "^#" $TEMP_DIR/$OUTBASE.$i.vcf > $TEMP_DIR/header.txt

    # merge the vcf region files
    if [[ "$ANNOTATE" == "true" ]]
    then
	FORK=""
	if [[ "$THREADS" -gt 1 ]]
	then
	    FORK="--fork $THREADS"
	fi

	if [[ $VERBOSE -eq 1 ]]
	then
	    echo -e "
            cat $TEMP_DIR/$OUTBASE.\"$chrom:\$start..\$end\".vcf | grep -v \"^#\" \\
                | $BEDTOOLS sort -i stdin | cat $TEMP_DIR/header.txt - \\
                | $VEP \\
                $FORK \\
                -o /dev/stdout \\
                --force_overwrite \\
                --offline \\
                --no_stats \\
                --cache \\
                --dir_cache $VEP_CACHE_DIR \\
                --species homo_sapiens \\
                --sift b \\
                --polyphen b \\
                --symbol \\
                --numbers \\
                --biotype \\
                --total_length \\
                --vcf \\
                --fields Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE \\
                | $BGZIP -c > $OUTPUT.vcf.gz"
	fi

	for i in `cat $WINDOWS | awk '{print $1":"$2".."$3}'`
	do
            # add "|| true" so it doesn't bail upon failure.
	    cat $TEMP_DIR/$OUTBASE."$i".vcf | grep -v "^#" || true
	done \
	    | $BEDTOOLS sort -i stdin | cat $TEMP_DIR/header.txt - \
	    | $VEP \
            $FORK \
            -o /dev/stdout \
	    --force_overwrite \
            --offline \
            --no_stats \
            --cache \
            --dir_cache $VEP_CACHE_DIR \
            --species homo_sapiens \
            --sift b \
            --polyphen b \
            --symbol \
            --numbers \
            --biotype \
            --total_length \
            --vcf \
            --fields Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE \
	    | $BGZIP -c > $OUTPUT.vcf.gz

    else
	if [[ $VERBOSE -eq 1 ]]
        then
            echo -e "
            cat $TEMP_DIR/$OUTBASE.\"$chrom:\$start..\$end\".vcf | grep -v \"^#\" \\
                | $BEDTOOLS sort -i stdin | cat $TEMP_DIR/header.txt - \\
                | $BGZIP -c > $OUTPUT.vcf.gz"
        fi

	for i in `cat $WINDOWS | awk '{print $1":"$2".."$3}'`
	do
            # if this fails then it bails out the script
	    cat $TEMP_DIR/$OUTBASE."$i".vcf | grep -v "^#" || true
	done \
	    | $BEDTOOLS sort -i stdin | cat $TEMP_DIR/header.txt - \
	    | $BGZIP -c > $OUTPUT.vcf.gz
    fi

    # index the vcf
    if [[ $VERBOSE -eq 1 ]]
    then
	echo -e "
	$TABIX -f -p vcf $OUTPUT.vcf.gz"
    fi
    $TABIX -f -p vcf $OUTPUT.vcf.gz

    # clean up
    if [[ "$KEEP" -eq 0 ]]
    then
	if [[ $VERBOSE -eq 1 ]]
	then
	    echo -e "
        rm -r $TEMP_DIR"
	fi
	rm -r $TEMP_DIR
    fi

    echo "Done"

    # exit cleanly
    exit 0

}

function somatic() {
    function somatic_usage() {
        echo "
usage:   speedseq somatic [options] <reference.fa> <normal.bam> <tumor.bam>

positional args:
         reference.fa
                  genome reference fasta file
         normal.bam
                  germline BAM file(s) (comma separated BAMs from multiple libraries).
                  Must have readgroup information, and the SM readgroup tag will
                  be the VCF column header
         tumor.bam
                  tumor BAM file(s) (comma separated BAMs for multiple libraries).
                    Must have readgroup information, and the SM readgroup tag will
                    be the VCF column header

options:
         -o STR   output prefix [tumor.bam]
         -w FILE  BED file of windowed genomic intervals
         -t INT   threads [1]
         -s       only output somatic variants
         -F FLOAT Require at least this fraction of observations supporting
                    an alternate allele within a single individual in order
                    to evaluate the position [0.05]
         -C INT   Require at least this count of observations supporting
                    an alternate allele within a single individual in order
                    to evaluate the position [2]
         -S FLOAT minimum somatic score (SSC) for PASS [18]
         -q FLOAT minimum QUAL score to output non-passing somatic variants [1e-5]
         -T DIR   temp directory [./output_prefix.XXXXXXXXXXXX]

         -A BOOL  annotate the vcf with VEP (true or false) (default: true)
         -K FILE  path to speedseq.config file (default: same directory as speedseq)
         -v       verbose
         -k       keep tempory files
         -h       show this message
"
    }

    # Check options passed in.
    if test -z "$3"
    then
        somatic_usage
        exit 1
    fi

    # set defaults
    SPEEDSEQ_DIR=`dirname $0`
    CONFIG="$SPEEDSEQ_DIR/speedseq.config"
    REF="${@:(-3):1}"
    OUTPUT=`basename "$TUMOR_BAM"`
    THREADS=1
    MIN_ALT_FRACTION=0.05
    MIN_ALT_COUNT=2
    TEMP_DIR=""
    ANNOTATE="true"
    MINQUAL=1e-5
    ONLY_SOMATIC=0
    SSC_THRES=18
    VERBOSE=0
    KEEP=0

    while getopts ":ho:w:t:F:C:T:A:q:S:svkK:" OPTION
    do
        case "${OPTION}" in
            h)
                somatic_usage
                exit 1
                ;;
            o)
                OUTPUT="$OPTARG"
                ;;
            w)
                WINDOWS="$OPTARG"
                ;;
            t)
                THREADS="$OPTARG"
                ;;
            F)
                MIN_ALT_FRACTION="$OPTARG"
                ;;
	    C)
		MIN_ALT_COUNT="$OPTARG"
		;;
            T)
                TEMP_DIR="$OPTARG"
                ;;
            A)
                ANNOTATE=`echo "$OPTARG" | tr [:upper:] [:lower:]`
                ;;
	    q)
		MINQUAL="$OPTARG"
		;;
	    S)
		SSC_THRES="$OPTARG"
		;;
	    s)
		ONLY_SOMATIC=1
		;;
            v)
                VERBOSE=1
                ;;
	    k)
		KEEP=1
		;;
	    K)
		CONFIG="$OPTARG"
		;;
        esac
    done

    NORMAL_BAM_STRING="${@:$((${OPTIND}+1)):1}"
    TUMOR_BAM_STRING="${@:$((${OPTIND}+2)):1}"
    TUMOR_BAM_LIST=(`echo $TUMOR_BAM_STRING | tr "," " "`)
    NORMAL_BAM_LIST=`echo $NORMAL_BAM_STRING | tr "," " "`

    if [[ -z $OUTPUT ]]
    then
	OUTPUT=`basename "${TUMOR_BAM_LIST[0]}"`
    fi
    OUTBASE=`basename "$OUTPUT"`
    OPTIND=0

    # Check the for the relevant binaries
    source_binaries $CONFIG

    if [[ ! -f "$FREEBAYES" ]]
    then
        somatic_usage
        echo -e "Error: freebayes executable not found. Please set path in $SPEEDSEQ_DIR/speedseq.config file\n"
        exit 1
    elif [[ ! -f "$BEDTOOLS" ]]
    then
        somatic_usage
        echo -e "Error: bedtools executable not found. Please set path in $SPEEDSEQ_DIR/speedseq.config file\n"
        exit 1
    elif [[ ! -f "$BGZIP" ]]
    then
        somatic_usage
        echo -e "Error: bgzip executable not found. Please set path in $SPEEDSEQ_DIR/speedseq.config file\n"
        exit 1
    elif [[ ! -f "$TABIX" ]]
    then
        somatic_usage
        echo -e "Error: tabix executable not found. Please set path in $SPEEDSEQ_DIR/speedseq.config file\n"
        exit 1
    elif [[ ! -f "$VEP" ]] && [[ "$ANNOTATE" == "true" ]]
    then
        somatic_usage
        echo -e "Error: VEP not found. Please set path in $SPEEDSEQ_DIR/speedseq.config file\n"
        exit 1
    elif [[ ! -d "$VEP_CACHE_DIR" ]] && [[ "$ANNOTATE" == "true" ]]
    then
        somatic_usage
        echo -e "Error: VEP cache directory not found. Please set path in $SPEEDSEQ_DIR/speedseq.config file\n"
        exit 1
    elif [[ ! -f "$PARALLEL" ]]
    then
        somatic_usage
        echo -e "Error: parallel executable not found. Please set path in $SPEEDSEQ_DIR/speedseq.config file\n"
        exit 1
    fi

    # Check that ANNOTATE is either true or false
    if [[ "$ANNOTATE" != "true" ]] && [[ "$ANNOTATE" != "false" ]]
    then
	somatic_usage
	echo -e "Error: -A must be either true or false\n"
	exit 1
    fi

    # Check that the ref and bam files exist
    if [[ -z "$REF" ]] || [[ ! -f "$REF" ]]
    then
        somatic_usage
        echo -e "Error: Reference file $REF not found.\n"
        exit 1
    fi
    for TEST_BAM in ${NORMAL_BAM_LIST[@]} ${TUMOR_BAM_LIST[@]}
    do
        if [[ ! -f $TEST_BAM ]]
        then
            somatic_usage
            echo -e "Error: BAM file $TEST_BAM not found.\n"
            exit 1
        fi
    done

    echo "Calling somatic variants..."
    # make temporary directory
    if [[ $VERBOSE -eq 1 ]]
    then
	echo "
        create temporary directory"
    fi
    if [[ -z $TEMP_DIR ]]
    then
	TEMP_DIR=`mktemp -d -p ./ ${OUTBASE}.XXXXXXXXXXXX`
    else
	mkdir -p $TEMP_DIR	
    fi

    # if no windows file, then make naive windows based on the chroms
    if [[ -z $WINDOWS ]]
    then
	if [[ $VERBOSE -eq 1 ]]
	then
	    echo -e "
        $SAMBAMBA view -H ${NORMAL_BAM_LIST[0]} | grep \"^@SQ\" | cut -f 2- | awk '{ gsub(\"^SN:\",\"\",\$1); gsub(\"^LN:\",\"\",\$2); print \$1\"\\\t0\\\t\"\$2; }' > $TEMP_DIR/windows.bed"
	fi
        $SAMBAMBA view -H ${NORMAL_BAM_LIST[0]} | grep "^@SQ" | cut -f 2- | awk '{ gsub("^SN:","",$1); gsub("^LN:","",$2); print $1"\t0\t"$2; }' > $TEMP_DIR/windows.bed
        WINDOWS="$TEMP_DIR/windows.bed"
    fi

    # write command to call variants on each of the windows in parallel
    NORMAL_BAM_FMT=`echo $NORMAL_BAM_STRING | tr "," " "`
    TUMOR_BAM_FMT=`echo $TUMOR_BAM_STRING | tr "," " "`

    if [[ $VERBOSE -eq 1 ]]
    then
	echo -e "
        $FREEBAYES -f $REF \\
            --pooled-discrete \\
            --min-repeat-entropy 1 \\
            --genotype-qualities \\
            --min-alternate-fraction ${MIN_ALT_FRACTION} \\
            --min-alternate-count ${MIN_ALT_COUNT} \\
            --region \$chrom:\$start..\$end \\
            $NORMAL_BAM_FMT $TUMOR_BAM_FMT \\
            | somatic_filter $MINQUAL $SSC_THRES $ONLY_SOMATIC \\
            > ${TEMP_DIR}/$OUTBASE.\$chrom:\$start..\$end.vcf"
    fi
    for i in `cat $WINDOWS | awk '{print $1":"$2".."$3}'`
    do
	# see the function somatic_filter for a more readable version of this awk command
	echo -e "$FREEBAYES -f $REF \
        --pooled-discrete \
        --genotype-qualities \
        --min-repeat-entropy 1 \
        --min-alternate-fraction ${MIN_ALT_FRACTION} \
        --min-alternate-count ${MIN_ALT_COUNT} \
        --region $i \
        $NORMAL_BAM_FMT $TUMOR_BAM_FMT \
        | awk -v ONLY_SOMATIC=\"$ONLY_SOMATIC\" -v MINQUAL=\"$MINQUAL\" -v SSC_THRES=\"$SSC_THRES\" 'BEGIN {NORMAL=10; TUMOR=11; GL_IDX=0;} { if (\$0~\"^#\") { print ; next; } if (! GL_IDX) { split(\$9,fmt,\":\") ; for (i=1;i<=length(fmt);++i) { if (fmt[i]==\"GL\") GL_IDX=i } } split(\$NORMAL,N,\":\"); split(N[GL_IDX],NGL,\",\"); split(\$TUMOR,T,\":\"); split(T[GL_IDX],TGL,\",\"); LOD_NORM=NGL[1]-NGL[2]; LOD_TUMOR_HET=TGL[2]-TGL[1]; LOD_TUMOR_HOM=TGL[3]-TGL[1]; if (LOD_TUMOR_HET > LOD_TUMOR_HOM) { LOD_TUMOR=LOD_TUMOR_HET } else { LOD_TUMOR=LOD_TUMOR_HOM } DQUAL=LOD_TUMOR+LOD_NORM; if (DQUAL>=SSC_THRES && \$NORMAL~\"^0/0\") { \$7=\"PASS\" ; \$8=\"SSC=\"DQUAL\";\"\$8 ; print } else if (!ONLY_SOMATIC && \$6>=MINQUAL && \$10~\"^0/0\" && ! match(\$11,\"^0/0\")) { \$8=\"SSC=\"DQUAL\";\"\$8 ; print } }' OFS=\"\t\" \
        > ${TEMP_DIR}/$OUTBASE.$i.vcf"
    done > $TEMP_DIR/var_command.txt

    # run the freebayes command in parallel
    if [[ $VERBOSE -eq 1 ]]
    then
	echo -e "
        cat $TEMP_DIR/var_command.txt | $PARALLEL -j $THREADS"
    fi
    cat $TEMP_DIR/var_command.txt | $PARALLEL -j $THREADS

    # make vcf header
    i=`head -n 1 $WINDOWS | awk '{print $1":"$2".."$3}'`
    if [[ $VERBOSE -eq 1 ]]
    then
	echo -e "
        grep \"^##\" $TEMP_DIR/$OUTBASE.$i.vcf \\
        | cat - <(echo '##INFO=<ID=SSC,Number=1,Type=Float,Description=\"Somatic score\">') <(grep \"^#CHROM\" $TEMP_DIR/$OUTBASE.$i.vcf) > $TEMP_DIR/header.txt"
    fi
    grep "^##" $TEMP_DIR/$OUTBASE.$i.vcf | cat - <(echo '##INFO=<ID=SSC,Number=1,Type=Float,Description="Somatic score">') <(grep "^#CHROM" $TEMP_DIR/$OUTBASE.$i.vcf) > $TEMP_DIR/header.txt

    # get the tumor and normal readgroups
    TUMOR_RG=`cat ${TEMP_DIR}/header.txt | tail -n 1 | cut -f 10`
    NORMAL_RG=`cat ${TEMP_DIR}/header.txt | tail -n 1 | cut -f 11`

    if [[ "$ANNOTATE" == "true" ]]
    then
        # merge the vcf region files, with VEP annotation
	FORK=""
	if [[ "$THREADS" -gt 1 ]]
	then
	    FORK="--fork $THREADS"
	fi

	if [[ $VERBOSE -eq 1 ]]
	then
	    echo -e "
            cat $TEMP_DIR/$OUTBASE.\"\$chrom:\$start..\$end\".vcf | grep -v \"^#\" \\
                | $BEDTOOLS sort -i stdin | cat $TEMP_DIR/header.txt - \\
                | $VEP \\
                $FORK \\
                -o /dev/stdout \\
                --force_overwrite \\
                --offline \\
                --no_stats \\
                --cache \\
                --dir_cache $VEP_CACHE_DIR \\
                --species homo_sapiens \\
                --sift b \\
                --polyphen b \\
                --symbol \\
                --numbers \\
                --biotype \\
                --total_length \\
                --vcf \\
                --fields Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE \\
                | $BGZIP -c > $OUTPUT.vcf.gz"
	fi    
	for i in `cat $WINDOWS | awk '{print $1":"$2".."$3}'`
	do
            # if this fails then it bails out the script, unless || true present
            cat $TEMP_DIR/$OUTBASE."$i".vcf | grep -v "^#" || true
	done \
            | $BEDTOOLS sort -i stdin | cat $TEMP_DIR/header.txt - \
            | $VEP \
            $FORK \
            -o /dev/stdout \
	    --force_overwrite \
            --offline \
            --no_stats \
            --cache \
            --dir_cache $VEP_CACHE_DIR \
            --species homo_sapiens \
            --sift b \
            --polyphen b \
            --symbol \
            --numbers \
            --biotype \
            --total_length \
            --vcf \
            --fields Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE \
            | $BGZIP -c > $OUTPUT.vcf.gz
    else
        # merge the vcf region files, without VEP annotation
	if [[ $VERBOSE -eq 1 ]]
	then
	    echo -e "
        cat $TEMP_DIR/$OUTBASE.\"\$chrom:\$start..\$end\".vcf | grep -v \"^#\" \\
            | $BEDTOOLS sort -i stdin | cat $TEMP_DIR/header.txt - \\
            | $BGZIP -c > $OUTPUT.vcf.gz"
	fi    
	for i in `cat $WINDOWS | awk '{print $1":"$2".."$3}'`
	do
            # if this fails then it bails out the script
            cat $TEMP_DIR/$OUTBASE."$i".vcf | grep -v "^#" || true
	done \
            | $BEDTOOLS sort -i stdin | cat $TEMP_DIR/header.txt - \
            | $BGZIP -c > $OUTPUT.vcf.gz
    fi

    # index the vcf
    if [[ $VERBOSE -eq 1 ]]
    then
	echo -e "
        $TABIX -f -p vcf $OUTPUT.vcf.gz"
    fi
    $TABIX -f -p vcf $OUTPUT.vcf.gz

    # produce PED file for GEMINI loading
    NORMAL_SAMPLE=`$SAMBAMBA view -H ${NORMAL_BAM_LIST[0]} | grep -m 1 "^@RG" | awk -v i=$i '{ for (j=1;j<=NF;++j) {if ($j~"^SM:") { gsub("^SM:","",$j); print $j } } }'`
    TUMOR_SAMPLE=`$SAMBAMBA view -H ${TUMOR_BAM_LIST[0]} | grep -m 1 "^@RG" | awk -v i=$i '{ for (j=1;j<=NF;++j) {if ($j~"^SM:") { gsub("^SM:","",$j); print $j } } }'`
    if [[ $VERBOSE -eq 1 ]]
    then
	echo "# Make PED file"
	echo "echo -e \"1\t$NORMAL_SAMPLE\tNone\tNone\t0\t1\n1\t$TUMOR_SAMPLE\tNone\tNone\t0\t2\" > $OUTPUT.ped"
    fi
    echo -e "1\t$NORMAL_SAMPLE\tNone\tNone\t0\t1\n1\t$TUMOR_SAMPLE\tNone\tNone\t0\t2" > $OUTPUT.ped

    # clean up
    if [[ "$KEEP" -eq 0 ]]
    then
	if [[ $VERBOSE -eq 1 ]]
	then
	    echo -e "
        rm -r $TEMP_DIR"
	fi
	rm -r $TEMP_DIR
    fi
	
    echo "Done"

    # exit cleanly
    exit 0
}

function sv() {
    function sv_usage() {
        echo "
usage:   speedseq sv [options]

sv options:
         -B FILE  full BAM file(s) (comma separated) (required)
         -S FILE  split reads BAM file(s) (comma separated, order as in -B) (required)
         -D FILE  discordant reads BAM files(s) (comma separated, order as in -B) (required)
         -R FILE  indexed reference genome fasta file (required)
         -o STR   output prefix [fullBam.bam]
         -t INT   threads [1] 
         -x FILE  BED file to exclude
         -g       genotype SV breakends with svtyper
         -d       calculate read-depth with CNVnator
         -A BOOL  annotate the vcf with VEP (true or false) (default: true)
         -m INT   minimum weight for a call [4]
         -r FLOAT trim threshold [0]
         -L INT   read length [auto]
         -T DIR   temp directory [./output_prefix.XXXXXXXXXXXX]
         -k       keep temporary files

         -s STR   lumpy split read parameters [auto]
                      bam_file:<splitreads.bam>,
                      back_distance:<10>,
                      min_mapping_threshold:<20>,
                      weight:<1>,
                      id:<11>,
                      min_clip:<20>

         -p STR   lumpy discordant read parameters [auto]
                      bam_file:<discreads.bam>,
                      histo_file:<auto>,
                      mean:<auto>,
                      stdev:<auto>,
                      read_length:<auto>,
                      min_non_overlap:<read_length>,
                      discordant_z:<5>,
                      back_distance:<10>,
                      min_mapping_threshold:<20>,
                      weight:<1>,
                      id:<10>

global options:
         -K FILE  path to speedseq.config file (default: same directory as speedseq)
         -v       verbose
         -h       show this message
"
    }

    # set defaults
    SPEEDSEQ_DIR=`dirname $0`
    CONFIG="$SPEEDSEQ_DIR/speedseq.config"
    THREADS=1
    ANNOTATE="true"
    MIN_WEIGHT=4
    TRIM_THRES=0
    EXCLUDE_BED=
    TEMP_DIR=""
    GENOTYPE=0
    READDEPTH=0
    VERBOSE=0
    KEEP=0
    OUTPUT=""

    while getopts ":hB:S:D:R:o:m:r:L:x:s:p:T:t:A:dgkvK:" OPTION
    do
	case "${OPTION}" in
            h)
                sv_usage
                exit 1
                ;;
	    B)
		FULL_BAM_STRING="$OPTARG"
		;;
	    S)
		SPL_BAM_STRING="$OPTARG"
		;;
	    D)
		DISC_BAM_STRING="$OPTARG"
		;;
	    R)
		REF="$OPTARG"
		;;
            o)
                OUTPUT="$OPTARG"
                ;;
            m)
                MIN_WEIGHT="$OPTARG"
                ;;
            r)
                TRIM_THRES="$OPTARG"
                ;;
	    L)
		READ_LENGTH="$OPTARG"
		;;
            x)
	        EXCLUDE_BED="$OPTARG"
		EXCLUDE_BED_FMT="-x $EXCLUDE_BED"
		;;  
	    s)
		LUMPY_SPL_STRING="-sr $OPTARG"
		;;
	    p)
		LUMPY_DISC_STRING="-pe $OPTARG"
		;;
            T)
                TEMP_DIR="$OPTARG"
                ;;
	    t)
		THREADS="$OPTARG"
		;;
	    A)
		ANNOTATE=`echo "$OPTARG" | tr [:upper:] [:lower:]`
		;;
	    d)
		READDEPTH=1
		;;
	    g)
		GENOTYPE=1
		;;
            v)
                VERBOSE=1
		;;
	    k)
		KEEP=1
		;;
	    K)
		CONFIG="$OPTARG"
		;;
	esac
    done

    # parse the BAM strings
    FULL_BAM_LIST=(`echo $FULL_BAM_STRING | tr "," " "`)
    SPL_BAM_LIST=(`echo $SPL_BAM_STRING | tr "," " "`)
    DISC_BAM_LIST=(`echo $DISC_BAM_STRING | tr "," " "`)

    OPTIND=0

    # Check the for the relevant binaries
    source_binaries $CONFIG

    if [[ -z "$LUMPY" ]]
    then
	sv_usage
        echo -e "Error: lumpy executable not found. Please set path in $SPEEDSEQ_DIR/speedseq.config file\n"
        exit 1
    elif [[ -z  "$PAIREND_DISTRO" ]]
    then
	sv_usage
        echo -e "Error: pairend_distro.py executable not found. Please set path in $SPEEDSEQ_DIR/speedseq.config file\n"
        exit 1
    elif [[ -z  "$SAMBAMBA" ]]
    then
	sv_usage
        echo -e "Error: sambamba executable not found. Please set path in $SPEEDSEQ_DIR/speedseq.config file\n"
        exit 1
    elif [[ -z  "$BEDPETOVCF" ]]
    then
	sv_usage
        echo -e "Error: bedpeToVcf executable not found. Please set path in $SPEEDSEQ_DIR/speedseq.config file\n"
        exit 1
    elif [[ ! -f "$VEP" ]] && [[ "$ANNOTATE" == "true" ]]
    then
	sv_usage
        echo -e "Error: VEP not found. Please set path in $SPEEDSEQ_DIR/speedseq.config file\n"
        exit 1
    elif [[ ! -d "$VEP_CACHE_DIR" ]] && [[ "$ANNOTATE" == "true" ]]
    then
	sv_usage
        echo -e "Error: VEP cache directory not found. Please set path in $SPEEDSEQ_DIR/speedseq.config file\n"
        exit 1
    fi

    # if genotyping requested, look for svtyper
    if [[ "$GENOTYPE" -eq 1 ]] && [[ -z "$SVTYPER" ]]
    then
	sv_usage
        echo -e "Error: svtyper executable not found. Please set path in $SPEEDSEQ_DIR/speedseq.config file\n"
        exit 1
    fi

    # if CNV read-depth requested, look for cnvnator executables
    if [[ "$READDEPTH" -eq 1 ]]
    then
	if [[ -z "$CNVNATOR_MULTI" ]]
	then
	    sv_usage
            echo -e "Error: cnvnator executable not found. Please set path in $SPEEDSEQ_DIR/speedseq.config file\n"
            exit 1
	elif [[ -z "$CNVNATOR_WRAPPER" ]]
	then
	    sv_usage
            echo -e "Error: cnvnator_wrapper.py executable not found. Please set path in $SPEEDSEQ_DIR/speedseq.config file\n"
            exit 1
	elif [[ -z "$ANNOTATE_RD" ]]
	then
	    sv_usage
            echo -e "Error: annotate_rd.py executable not found. Please set path in $SPEEDSEQ_DIR/speedseq.config file\n"
            exit 1
	fi
    fi

    # check for required python modules (pysam, numpy, scipy, etc)
    check_python_modules $PYTHON27

    # Check that the required files exist
    if [[ ! -f $REF ]]
    then
	sv_usage
	echo -e "Error: reference fasta file $REF not found\n"
	exit 1
    fi
    if [[ ! -f $REF.fai && ! -f $(echo ${REF%*.*}).fai ]]
    then
	sv_usage
	echo -e "Error: reference fasta file $REF not indexed\n"
	exit 1
    fi
    if [[ ${#FULL_BAM_LIST[@]} -eq 0 ]] || [[ ${#SPL_BAM_LIST[@]} -eq 0 ]] || [[ ${#DISC_BAM_LIST[@]} -eq 0 ]]
    then
	sv_usage
	echo -e "Error: -B, -S, and -D are required\n"
	exit 1
    fi
    if [[ ! ${#FULL_BAM_LIST[@]} -eq ${#SPL_BAM_LIST[@]} ]] || [[ ! ${#FULL_BAM_LIST[@]} -eq ${#DISC_BAM_LIST[@]} ]]
    then
	sv_usage
	echo -e "Error: -B, -S, and -D arguments have different lengths\n"
	exit 1
    fi
    for TEST_BAM in ${FULL_BAM_LIST[@]} ${SPL_BAM_LIST[@]} ${DISC_BAM_LIST[@]}
    do
        if [[ ! -f $TEST_BAM ]]
        then
            sv_usage
            echo -e "Error: BAM file $TEST_BAM not found.\n"
            exit 1
        fi
    done

    # Check that ANNOTATE is either true or false
    if [[ "$ANNOTATE" != "true" ]] && [[ "$ANNOTATE" != "false" ]]
    then
	sv_usage
	echo -e "Error: -A must be either true or false\n"
	exit 1
    fi

    # default OUTPUT if not provided
    if test -z "$OUTPUT"
    then
	OUTPUT=`basename "${FULL_BAM_LIST[0]}"`
    fi
    OUTBASE=`basename "$OUTPUT"`

    # make temporary directory
    if [[ $VERBOSE -eq 1 ]]
    then
	echo "
        create temporary directory"
    fi
    if [[ -z $TEMP_DIR ]]
    then
	TEMP_DIR=`mktemp -d -p ./ ${OUTBASE}.XXXXXXXXXXXX`
    else
	mkdir -p $TEMP_DIR	
    fi

    # create the LUMPY config file for generating the VCF
    for i in $( seq 0 $(( ${#FULL_BAM_LIST[@]}-1 )) )
    do
        FULL_BAM=${FULL_BAM_LIST[$i]}
        SPL_BAM=${SPL_BAM_LIST[$i]}
        DISC_BAM=${DISC_BAM_LIST[$i]}

	FULL_SAMPLE=`$SAMBAMBA view -H $FULL_BAM | grep -m 1 "^@RG" | awk -v i=$i '{ for (j=1;j<=NF;++j) {if ($j~"^SM:") { gsub("^SM:","",$j); print $j } } }'`
	SPL_SAMPLE=`$SAMBAMBA view -H $SPL_BAM | grep -m 1 "^@RG" | awk -v i=$i '{ for (j=1;j<=NF;++j) {if ($j~"^SM:") { gsub("^SM:","",$j); print $j } } }'`
	DISC_SAMPLE=`$SAMBAMBA view -H $DISC_BAM | grep -m 1 "^@RG" | awk -v i=$i '{ for (j=1;j<=NF;++j) {if ($j~"^SM:") { gsub("^SM:","",$j); print $j } } }'`

	echo -e "$DISC_SAMPLE\t$(($i+1))0\tPE\t$DISC_BAM"
	echo -e "$SPL_SAMPLE\t$(($i+1))1\tSR\t$SPL_BAM"
	echo -e "$FULL_SAMPLE\t$(($i+1))2\tRD\t$FULL_BAM"
	
    done > $TEMP_DIR/$OUTBASE.sample.config

    # if LUMPY_SPL_STRING not provided, make one
    if test -z "$LUMPY_SPL_STRING"
    then
	for i in $( seq 0 $(( ${#FULL_BAM_LIST[@]}-1 )) )
	do
            FULL_BAM=${FULL_BAM_LIST[$i]}
            SPL_BAM=${SPL_BAM_LIST[$i]}
            DISC_BAM=${DISC_BAM_LIST[$i]}
	    
            FULL_BASE=`basename "$FULL_BAM"`
            SPL_BASE=`basename "$SPL_BAM"`
            DISC_BASE=`basename "$DISC_BAM"`

	    LUMPY_SPL_STRING="$LUMPY_SPL_STRING -sr bam_file:${SPL_BAM},back_distance:10,min_mapping_threshold:20,weight:1,id:$(($i+1))1,min_clip:20"
	done
    fi

    # if LUMPY_DISC_STRING not provided, make one
    if test -z "$LUMPY_DISC_STRING"
    then
	echo -n "Calculating alignment stats... "
	for i in $( seq 0 $(( ${#FULL_BAM_LIST[@]}-1 )) )
	do
            FULL_BAM=${FULL_BAM_LIST[$i]}
            SPL_BAM=${SPL_BAM_LIST[$i]}
            DISC_BAM=${DISC_BAM_LIST[$i]}
	    
            FULL_BASE=`basename "$FULL_BAM"`
            SPL_BASE=`basename "$SPL_BAM"`
            DISC_BASE=`basename "$DISC_BAM"`
	    
            # calc readlength if not provided
            if test -z "$READ_LENGTH"
            then
		READ_LENGTH=`$SAMBAMBA view $FULL_BAM | head -n 10000 | awk 'BEGIN { MAX_LEN=0 } { LEN=length($10); if (LEN>MAX_LEN) MAX_LEN=LEN } END { print MAX_LEN }'`
            fi

	    if [[ "VERBOSE" -eq 1 ]]
	    then
		echo "
# calculate stats
$SAMBAMBA view -F \"paired and mate_is_reverse_strand and not (unmapped or mate_is_unmapped or reverse_strand or secondary_alignment or duplicate)\" $FULL_BAM \\
    | head -n 10000000 | tail -n 1000000 \\
    | $PAIREND_DISTRO -r $READ_LENGTH -X 4 -N 1000000 -o ${TEMP_DIR}/$FULL_BASE.x4.histo \\
    > ${TEMP_DIR}/$FULL_BASE.insert.stats
	    "
	    fi

            # calculate stats
            $SAMBAMBA view -F "paired and mate_is_reverse_strand and not (unmapped or mate_is_unmapped or reverse_strand or secondary_alignment or duplicate)" $FULL_BAM \
		| head -n 10000000 | tail -n 1000000 \
		| $PAIREND_DISTRO -r $READ_LENGTH -X 4 -N 1000000 -o ${TEMP_DIR}/$FULL_BASE.x4.histo \
		> ${TEMP_DIR}/$FULL_BASE.insert.stats
	    
            MEAN=`cat ${TEMP_DIR}/$FULL_BASE.insert.stats | tr '\t' '\n' | grep "^mean" | sed 's/mean\://g'`
            STDEV=`cat ${TEMP_DIR}/$FULL_BASE.insert.stats | tr '\t' '\n' | grep "^stdev" | sed 's/stdev\://g'`

	    LUMPY_DISC_STRING="$LUMPY_DISC_STRING -pe bam_file:${DISC_BAM},histo_file:${TEMP_DIR}/${FULL_BASE}.x4.histo,mean:${MEAN},stdev:${STDEV},read_length:${READ_LENGTH},min_non_overlap:${READ_LENGTH},discordant_z:5,back_distance:10,weight:1,id:$(($i+1))0,min_mapping_threshold:20"
	done
	echo "done"
    fi

    echo "Running LUMPY... "
    if [[ "$VERBOSE" -eq 1 ]]
    then
	echo "
# call lumpy
$LUMPY -t ${TEMP_DIR}/${OUTBASE} -mw $MIN_WEIGHT -tt $TRIM_THRES $EXCLUDE_BED_FMT $LUMPY_DISC_STRING $LUMPY_SPL_STRING > $TEMP_DIR/$OUTBASE.sv.lumpy
        "
    fi
    # call lumpy
    $LUMPY -t ${TEMP_DIR}/${OUTBASE} -mw $MIN_WEIGHT -tt $TRIM_THRES \
    	$EXCLUDE_BED_FMT \
    	$LUMPY_DISC_STRING \
    	$LUMPY_SPL_STRING \
    	> $TEMP_DIR/$OUTBASE.sv.lumpy

    if [[ "$VERBOSE" -eq 1 ]]
    then
	echo -e "
# convert to VCF, ensure that each variant has at least one sample with MIN_WEIGHT support
cat $TEMP_DIR/$OUTBASE.sv.lumpy \\
    | $PYTHON27 $BEDPETOVCF -t LUMPY -c $TEMP_DIR/$OUTBASE.sample.config -f $REF \\
    | awk -v MIN_WEIGHT=$MIN_WEIGHT 'BEGIN {S_IDX=0} {
        if (\$0~\"^#\") print
        else {
            if (S_IDX==0) {
                split(\$9,fmt,\":\");
                for (i=1;i<=length(fmt);++i) {
                    if (fmt[i]==\"SUP\") S_IDX=i
                }
            }
            if (S_IDX!=0) {
                PASS=0
                for (i=10;i<=NF;++i) {
                    split(\$i,a,\":\");
                    if (a[S_IDX]>=MIN_WEIGHT) PASS=1;
                }
            }
        }
        if (PASS) print
    }' \\
    > $TEMP_DIR/$OUTBASE.sv.vcf
        "
    fi

    # convert to VCF, ensure that each variant has at least one sample with MIN_WEIGHT support
    cat $TEMP_DIR/$OUTBASE.sv.lumpy \
	| $PYTHON27 $BEDPETOVCF -t LUMPY -c $TEMP_DIR/$OUTBASE.sample.config -f $REF \
	| awk -v MIN_WEIGHT=$MIN_WEIGHT 'BEGIN {S_IDX=0} {
            if ($0~"^#") print
            else {
                if (S_IDX==0) {
                    split($9,fmt,":");
                    for (i=1;i<=length(fmt);++i) {
                        if (fmt[i]=="SUP") S_IDX=i
                    }
                }
                if (S_IDX!=0) {
                    PASS=0
                    for (i=10;i<=NF;++i) {
                        split($i,a,":");
                        if (a[S_IDX]>=MIN_WEIGHT) PASS=1;
                    }
                }
            }
            if (PASS) print
        }' \
	> $TEMP_DIR/$OUTBASE.sv.vcf

    if [[ "$VERBOSE" -eq 1 ]]
    then
	echo -e "
# make header
cat $TEMP_DIR/$OUTBASE.sv.vcf | awk '$0!~\"^#\" {exit}; 1' > $TEMP_DIR/header.txt
        "
    fi
    # make header
    cat $TEMP_DIR/$OUTBASE.sv.vcf | awk '$0!~"^#" {exit}; 1' > $TEMP_DIR/header.txt

    # genotype
    # WARNING: svtyper only works when there is one BAM per sample, not multi-library samples
    if [[ GENOTYPE -eq 1 ]]
    then
	for i in $( seq 0 $(( ${#FULL_BAM_LIST[@]}-1 )) )
	do
	    FULL_BAM=${FULL_BAM_LIST[$i]}
	    SPL_BAM=${SPL_BAM_LIST[$i]}
	    FULL_BASE=`basename "$FULL_BAM"`

	    MEAN=`cat ${TEMP_DIR}/$FULL_BASE.insert.stats | tr '\t' '\n' | grep "^mean" | sed 's/mean\://g'`
	    STDEV=`cat ${TEMP_DIR}/$FULL_BASE.insert.stats | tr '\t' '\n' | grep "^stdev" | sed 's/stdev\://g'`
	    READ_LENGTH=`$SAMBAMBA view $FULL_BAM | head -n 10000 | awk 'BEGIN { MAX_LEN=0 } { LEN=length($10); if (LEN>MAX_LEN) MAX_LEN=LEN } END { print MAX_LEN }'`

	    if [[ "$VERBOSE" -eq 1 ]]
	    then
		echo "# genotype structural variants"
		echo -e "$PYTHON27 $SVTYPER -v $TEMP_DIR/$OUTBASE.sv.vcf -B $FULL_BAM -S $SPL_BAM -m $MEAN -sd $STDEV -rl $READ_LENGTH > \$vcf.gt ; mv \$vcf.gt \$vcf"
	    fi

	    # genotype structural variants
	    $PYTHON27 $SVTYPER -v $TEMP_DIR/$OUTBASE.sv.vcf -B $FULL_BAM -S $SPL_BAM -m $MEAN -sd $STDEV -rl $READ_LENGTH > $TEMP_DIR/$OUTBASE.sv.gt.vcf ; mv $TEMP_DIR/$OUTBASE.sv.gt.vcf $TEMP_DIR/$OUTBASE.sv.vcf
	done
    fi

    # run cnvnator
    # WARNING: CNVNATOR currently only works when there is one BAM per sample, not multi-library samples
    if [[ "$READDEPTH" -eq 1 ]]
    then
	echo "Calculating read depth"
	WINDOW_SIZE=100

	for i in $( seq 0 $(( ${#FULL_BAM_LIST[@]}-1 )) )
	do
	    FULL_BAM=${FULL_BAM_LIST[$i]}
	    SAMPLE=`cat $TEMP_DIR/$OUTBASE.sample.config | awk -v FULL_BAM=$FULL_BAM '{ if ($4==FULL_BAM) { print $1; exit } }'`

	    if [[ "$VERBOSE" -eq 1 ]]
	    then
		echo "
    # run cnvnator-multi
    $PYTHON27 $CNVNATOR_WRAPPER --cnvnator $CNVNATOR_MULTI -T $TEMP_DIR/cnvnator-temp -t $THREADS -w 100 -b ${FULL_BAM} -o $TEMP_DIR/$OUTBASE.rd -c $CNVNATOR_CHROMS_DIR -g GRCh37
		"
	    fi
	    $PYTHON27 $CNVNATOR_WRAPPER --cnvnator $CNVNATOR_MULTI -T $TEMP_DIR/cnvnator-temp -t $THREADS -w 100 -b ${FULL_BAM} -o $TEMP_DIR/$OUTBASE.readdepth -c $CNVNATOR_CHROMS_DIR -g GRCh37

	    # Calculate read-depth of LUMPY calls
	    if [[ "$VERBOSE" -eq 1 ]]
	    then
		echo "
    # Calculate read-depth of LUMPY calls
    $PYTHON27 $ANNOTATE_RD --cnvnator $CNVNATOR_MULTI -s $SAMPLE -w $WINDOW_SIZE -r $TEMP_DIR/cnvnator-temp/${FULL_BAM}.hist.root -v $TEMP_DIR/$OUTBASE.sv.vcf > $TEMP_DIR/$OUTBASE.sv.rd.vcf
		"
	    fi
	    $PYTHON27 $ANNOTATE_RD --cnvnator $CNVNATOR_MULTI -s $SAMPLE -w $WINDOW_SIZE -r $TEMP_DIR/cnvnator-temp/${FULL_BAM}.hist.root -v $TEMP_DIR/$OUTBASE.sv.vcf > $TEMP_DIR/$OUTBASE.rd.sv.vcf

	    if [[ "$VERBOSE" -eq 1 ]]
	    then
		echo "mv $TEMP_DIR/$OUTBASE.rd.sv.vcf $TEMP_DIR/$OUTBASE.rd.sv.vcf"
	    fi
	    mv $TEMP_DIR/$OUTBASE.rd.sv.vcf $TEMP_DIR/$OUTBASE.sv.vcf
	done
    fi

    # Annotate structural variants with VEP
    if [[ "$ANNOTATE" == "true" ]]
    then
	FORK=""
	if [[ "$THREADS" -gt 1 ]]
	then
	    FORK="--fork $THREADS"
	fi

	if [[ "$VERBOSE" -eq 1 ]]
	then
	    echo -e "
# Annotate structural variants with VEP
cat  $TEMP_DIR/$OUTBASE.sv.vcf \\
    | awk '\$0~\"^#\" || (\$1<=22 || \$1==\"X\" || \$1==\"Y\")' \\
    | $VEP \\
    $FORK \\
    -o /dev/stdout \\
    --force_overwrite \\
    --format vcf \\
    --offline \\
    --no_stats \\
    --cache \\
    --dir_cache $VEP_CACHE_DIR \\
    --species homo_sapiens \\
    --sift b \\
    --polyphen b \\
    --symbol \\
    --numbers \\
    --biotype \\
    --total_length \\
    --vcf \\
    --fields Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE \\
    | cat - <(cat $TEMP_DIR/$OUTBASE.sv.vcf | grep -v \"^#\" | awk '(\$1>22 && \$1!=\"X\" && \$1!=\"Y\")') \\
    > $TEMP_DIR/$OUTBASE.vep.sv.vcf
    "
	fi

	cat  $TEMP_DIR/$OUTBASE.sv.vcf \
	    | awk '$0~"^#" || (($1<=22 || $1=="X" || $1=="Y"))' \
	    | $VEP \
            $FORK \
            -o /dev/stdout \
            --force_overwrite \
	    --format vcf \
            --offline \
            --no_stats \
            --cache \
            --dir_cache $VEP_CACHE_DIR \
            --species homo_sapiens \
            --sift b \
            --polyphen b \
            --symbol \
            --numbers \
            --biotype \
            --total_length \
            --vcf \
            --fields Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE \
	    | cat - <(cat $TEMP_DIR/$OUTBASE.sv.vcf | grep -v "^#" | awk '($1>22 && $1!="X" && $1!="Y")') \
	    > $TEMP_DIR/$OUTBASE.vep.sv.vcf

	mv $TEMP_DIR/$OUTBASE.vep.sv.vcf $TEMP_DIR/$OUTBASE.sv.vcf
    fi


    # copy raw files to output directory
    cp $TEMP_DIR/$OUTBASE.sv.lumpy $OUTPUT.sv.lumpy
    cp $TEMP_DIR/$OUTBASE.sample.config $OUTPUT.sample.config

    # write output vcf file
    $BEDTOOLS sort -header -i $TEMP_DIR/$OUTBASE.sv.vcf \
	| $BGZIP -c \
	> $OUTPUT.sv.vcf.gz
    $TABIX -p vcf $OUTPUT.sv.vcf.gz

    # write output bedpe file
    $LUMPYTOBEDPE -b $TEMP_DIR/$OUTBASE.sv.lumpy -c $TEMP_DIR/$OUTBASE.sample.config \
	| awk -v MIN_WEIGHT=${MIN_WEIGHT} '{ PASS=0 ; for (i=16;i<=NF;++i) { if ($i>=MIN_WEIGHT) { PASS=1 ; print ; break } } }' \
	> $OUTPUT.sv.bedpe

    # clean up
    if [[ "$KEEP" -eq 0 ]]
    then
	rm -r ${TEMP_DIR}
    fi

    echo "done"

    # exit cleanly
    exit 0
}

function realign() {
    function realign_usage() {
	echo "
usage:   speedseq realign [options] <reference.fa> <in.bam>

positional args:
         reference.fa
                  fasta file (indexed with bwa)
         in.bam   re-align from a coordinate sorted BAM file

alignment options:
         -o STR   output prefix [in.realign]
         -I FLOAT[,FLOAT[,INT[,INT]]]
                  specify the mean, standard deviation (10% of the mean if absent), max
                    (4 sigma from the mean if absent) and min of the insert size distribution.
                    FR orientation only. [inferred]
         -g STR   re-align these read groups from the coordinate sorted BAM file (comma separated)
                    (default: re-align all readgroups)
         -n       rename reads for smaller file size
         -e       merge BAM files from different libraries
         -t INT   threads [1]
         -T DIR   temp directory [./output_prefix.XXXXXXXXXXXX]

samblaster options:
         -i       include duplicates in splitters and discordants
         -c INT   maximum number of split alignments for a read to be included in splitter file [2]
         -m INT   minimum non-overlapping base pairs between two alignments for a read to be included in splitter file [20]

sambamba options:
         -M       amount of memory in GB to be used for sorting [20]

global options:
         -K FILE  path to speedseq.config file (default: same directory as speedseq)
         -v       verbose
         -h       show this message
"
    }

    # Check options passed in.
    if test -z "$2"
    then
	realign_usage
	exit 1
    fi

    # set defaults
    SPEEDSEQ_DIR=`dirname $0`
    CONFIG="$SPEEDSEQ_DIR/speedseq.config"
    INTERLEAVED=0
    RG_FMT=""
    OUTPUT=""
    INCLUDE_DUPS="--excludeDups"
    MAX_SPLIT_COUNT=2
    MIN_NON_OVERLAP=20
    THREADS=1
    TEMP_DIR=""
    VERBOSE=0
    INS_DIST=""
    REALIGN_RG_LIST=""
    RENAME=""
    MERGE=0
    SORT_MEM=20 # amount of memory for sorting, in gigabytes

    while getopts ":hw:o:R:pic:m:M:t:T:g:I:nevK:" OPTION
    do
	case "${OPTION}" in
	    h)
		realign_usage
		exit 1
		;;
	    R)
		RG="$OPTARG"
		RG_FMT="-R '$OPTARG'"
		;;
	    p)
		INTERLEAVED=1
		;;
	    o)
		OUTPUT="$OPTARG"
		;;
	    i)
		INCLUDE_DUPS=""
		;;
	    c)
		MAX_SPLIT_COUNT="$OPTARG"
		;;
	    m)
		MIN_NON_OVERLAP="$OPTARG"
		;;
	    M)
		SORT_MEM="$OPTARG"
		;;
	    n)
		RENAME="-n"
		;;
	    t)
		THREADS="$OPTARG"
		;;
	    T)
		TEMP_DIR="$OPTARG"
		;;
	    g)
		REALIGN_RG_LIST="$OPTARG"
		;;
	    I)
		INS_DIST="-I $OPTARG"
		;;
	    v)
		VERBOSE=1
		;;
	    e)
		MERGE=1
		;;
	    K)
		CONFIG="$OPTARG"
		;;
	esac
    done

    REF="${@:${OPTIND}:1}"
    IN_BAM="${@:$((${OPTIND}+1)):1}"
    if [[ -z "$OUTPUT" ]]
    then
	OUTPUT=`basename "$IN_BAM" ".bam"`".realign"
    fi

    # Check that the ref and bam files exist
    if [[ -z "$REF" ]] || [[ ! -f "$REF" ]]
    then
	realign_usage
	echo -e "Error: Reference file $REF not found.\n"
	exit 1
    elif [[ -z "$IN_BAM" ]] || [[ ! -e "$IN_BAM" ]]
    then
	realign_usage
	echo -e "Error: BAM file $IN_BAM not found.\n"
	exit 1
    fi

    # Check the for the relevant binaries
    source_binaries $CONFIG
    if [[ -z "$BWA" ]]
    then
	realign_usage
        echo -e "Error: bwa executable not found. Please set path in $SPEEDSEQ_DIR/speedseq.config file\n"
        exit 1
    elif [[ -z  "$SAMBLASTER" ]]
    then
	realign_usage
        echo -e "Error: samblaster executable not found. Please set path in $SPEEDSEQ_DIR/speedseq.config file\n"
        exit 1
    elif [[ -z "$SAMBAMBA" ]]
    then
	realign_usage
        echo -e "Error: sambamba executable not found. Please set path in $SPEEDSEQ_DIR/speedseq.config file\n"
        exit 1
    elif [[ -z "$MBUFFER" ]]
    then
	realign_usage
	echo -e "Error: mbuffer executable not found. Please set path in $SPEEDSEQ_DIR/speedseq.config file\n"
	exit 1
    elif [[ -z "$BAMTOFASTQ" ]]
    then
	realign_usage
	echo -e "Error: bamtofastq.py executable not found. Please set path in $SPEEDSEQ_DIR/speedseq.config file\n"
	exit 1
    elif [[ -z "$BAMHEADRG" ]]
    then
	realign_usage
	echo -e "Error: bamheadrg.py executable not found. Please set path in $SPEEDSEQ_DIR/speedseq.config file\n"
	exit 1
    elif [[ -z "$PARALLEL" ]]
    then
        realign_usage
	echo parallel is $PARALLEL
        echo -e "Error: parallel executable not found. Please set path in $SPEEDSEQ_DIR/speedseq.config file\n"
        exit 1
    fi

    # check for required python modules (pysam, numpy, scipy, etc)
    check_python_modules $PYTHON27

    # parse the libraries in the BAM header to extract readgroups from the same library
    if [[ -z $REALIGN_RG_LIST ]]
    then
	REALIGN_RG_LIST=(`$SAMBAMBA view -H $IN_BAM \
	    | grep "^@RG" \
	    | awk -v i=$i \
		'{
		    for (j=1;j<=NF;++j) {
			if ($j~"^ID") {
			    id=$j
			    gsub("^ID:","",id)
			}
			else if ($j~"^LB:") {
			    lib=$j
			    gsub("^LB:","",lib);
			}
		    }
		    lib_list[lib][length(lib_list[lib])+1] = id
		}
		END {
		    for (i in lib_list) {
			RG_STRING=lib_list[i][1]
			for (j=2;j<=length(lib_list[i]);++j) {
			    RG_STRING=RG_STRING","lib_list[i][j]
			}
			print RG_STRING
		    }
		}'`)
    fi

    # set the output name
    OUTBASE=`basename "$OUTPUT"`

    echo "Aligning..."
    # create temp directory if not specified by command argument
    if [[ -z $TEMP_DIR ]]
    then
	TEMP_DIR=`mktemp -d -p ./ ${OUTBASE}.XXXXXXXXXXXX`
    fi

    if [[ $VERBOSE -eq 1 ]]
    then
	echo "
        mkdir -p $TEMP_DIR/full $TEMP_DIR/spl $TEMP_DIR/disc
        mkfifo $TEMP_DIR/spl_pipe $TEMP_DIR/disc_pipe"
    fi

    if [[ $VERBOSE -eq 1 ]]
    then
	echo "
        mkdir -p $TEMP_DIR $TEMP_DIR/full $TEMP_DIR/spl $TEMP_DIR/disc
        mkdir -p $TEMP_DIR/full $TEMP_DIR/spl $TEMP_DIR/disc
        mkfifo $TEMP_DIR/spl_pipe $TEMP_DIR/disc_pipe $TEMP_DIR/fq_pipe"
    fi

    # create temp files and pipes
    mkdir -p $TEMP_DIR $TEMP_DIR/full $TEMP_DIR/spl $TEMP_DIR/disc
    if [[ ! -e $TEMP_DIR/spl_pipe ]]
    then
	mkfifo $TEMP_DIR/spl_pipe
    fi
    if [[ ! -e $TEMP_DIR/disc_pipe ]]
    then
	mkfifo $TEMP_DIR/disc_pipe
    fi
    if [[ ! -e $TEMP_DIR/fq_pipe ]]
    then
	mkfifo $TEMP_DIR/fq_pipe
    fi
    FQ="$TEMP_DIR/fq_pipe"

    for i in $( seq 0 $(( ${#REALIGN_RG_LIST[@]}-1 )) )
    do
	echo ${REALIGN_RG_LIST[i]}
	REALIGN_RG=${REALIGN_RG_LIST[i]}
	REALIGN_RG_FMT="-r $REALIGN_RG"

	# create the virtual file for fastq pipe
	if [[ $VERBOSE -eq 1 ]]
	then
	    echo "
	    $PYTHON27 $BAMTOFASTQ $REALIGN_RG_FMT $RENAME -i $IN_BAM | $MBUFFER -q -m 1G > $TEMP_DIR/fq_pipe &"
	fi
	$PYTHON27 $BAMTOFASTQ $REALIGN_RG_FMT $RENAME -i $IN_BAM | $MBUFFER -q -m 1G > $TEMP_DIR/fq_pipe &

	# alignment command
	if [[ $VERBOSE -eq 1 ]]
	then
	    echo -e "
	$BWA mem -t $THREADS -C -M -p $INS_DIST $RG_FMT $REF $FQ | \\
	    $PYTHON27 $BAMHEADRG -b $IN_BAM $REALIGN_RG_FMT | \\
	    $SAMBLASTER $INCLUDE_DUPS --addMateTags --maxSplitCount $MAX_SPLIT_COUNT --minNonOverlap $MIN_NON_OVERLAP --splitterFile $TEMP_DIR/spl_pipe --discordantFile $TEMP_DIR/disc_pipe | \\
	    $SAMBAMBA view -S -f bam -l 0 /dev/stdin | \\
	    $SAMBAMBA sort -t $THREADS -m $((${SORT_MEM}-2))G --tmpdir=$TEMP_DIR/full -o $OUTPUT.$(($i+1)).bam /dev/stdin

	$SAMBAMBA view -S -f bam -l 0 $TEMP_DIR/spl_pipe | \\
	    $SAMBAMBA sort -t 4 -m 1G --tmpdir=$TEMP_DIR/spl -o $OUTPUT.$(($i+1)).splitters.bam /dev/stdin
	$SAMBAMBA view -S -f bam $TEMP_DIR/disc_pipe | \\
	    $SAMBAMBA sort -t 4 -m 1G --tmpdir=$TEMP_DIR/disc -o $OUTPUT.$(($i+1)).discordants.bam /dev/stdin"
	fi

	echo -e "
	$BWA mem -t $THREADS -C -M -p $INS_DIST $RG_FMT $REF $FQ | \
	    $PYTHON27 $BAMHEADRG -b <($SAMBAMBA view -f bam -l 0 $IN_BAM) $REALIGN_RG_FMT | \
	    $SAMBLASTER $INCLUDE_DUPS --addMateTags --maxSplitCount $MAX_SPLIT_COUNT --minNonOverlap $MIN_NON_OVERLAP --splitterFile $TEMP_DIR/spl_pipe --discordantFile $TEMP_DIR/disc_pipe | \
	    $SAMBAMBA view -S -f bam -l 0 /dev/stdin | \
	    $SAMBAMBA sort -t $THREADS -m $((${SORT_MEM}-2))G --tmpdir=$TEMP_DIR/full -o $OUTPUT.$(($i+1)).bam /dev/stdin

	$SAMBAMBA view -S -f bam -l 0 $TEMP_DIR/spl_pipe | \
	    $SAMBAMBA sort -t 4 -m 1G --tmpdir=$TEMP_DIR/spl -o $OUTPUT.$(($i+1)).splitters.bam /dev/stdin
	$SAMBAMBA view -S -f bam $TEMP_DIR/disc_pipe | \
	    $SAMBAMBA sort -t 4 -m 1G --tmpdir=$TEMP_DIR/disc -o $OUTPUT.$(($i+1)).discordants.bam /dev/stdin
	" | $PARALLEL -j 3

	# index the files unless they will be merged
	if [[ $MERGE -eq 0 ]] && [[ ${#REALIGN_RG_LIST[@]} -gt 1 ]]
	then
	    if [[ $VERBOSE -eq 1 ]]
	    then
		echo -e "
		$SAMBAMBA index $OUTPUT.$(($i+1)).bam
		$SAMBAMBA index $OUTPUT.$(($i+1)).discordants.bam
		$SAMBAMBA index $OUTPUT.$(($i+1)).splitters.bam"
	    fi

	    echo "
	    $SAMBAMBA index $OUTPUT.$(($i+1)).bam
	    $SAMBAMBA index $OUTPUT.$(($i+1)).discordants.bam
	    $SAMBAMBA index $OUTPUT.$(($i+1)).splitters.bam
	    " | $PARALLEL -j 3
	fi
    done

    # if only 1 library, then rename files
    if [[ ${#REALIGN_RG_LIST[@]} -eq 1 ]]
    then
	mv $OUTPUT.1.bam $OUTPUT.bam
	mv $OUTPUT.1.discordants.bam $OUTPUT.discordants.bam
	mv $OUTPUT.1.splitters.bam $OUTPUT.splitters.bam

	echo "
	$SAMBAMBA index $OUTPUT.bam
	$SAMBAMBA index $OUTPUT.discordants.bam
	$SAMBAMBA index $OUTPUT.splitters.bam
	" | $PARALLEL -j 3

    fi

    # merge and index the files, if requested
    if [[ $MERGE -eq 1 ]] && [[ ${#REALIGN_RG_LIST[@]} -gt 1 ]]
    then
	MERGE_FULL=""
	MERGE_DISCORDANTS=""
	MERGE_SPLITTERS=""
	for i in $( seq 0 $(( ${#REALIGN_RG_LIST[@]}-1 )) )
	do
	    MERGE_FULL="$MERGE_FULL $OUTPUT.$(($i+1)).bam"
	    MERGE_DISCORDANTS="$MERGE_DISCORDANTS $OUTPUT.$(($i+1)).discordants.bam"
	    MERGE_SPLITTERS="$MERGE_SPLITTERS $OUTPUT.$(($i+1)).splitters.bam"
	done

	if [[ $VERBOSE -eq 1 ]]
	then
	    echo "
            $SAMBAMBA merge -t $THREADS $OUTPUT.bam $MERGE_FULL
            $SAMBAMBA merge -t $THREADS $OUTPUT.discordants.bam $MERGE_DISCORDANTS
            $SAMBAMBA merge -t $THREADS $OUTPUT.splitters.bam $MERGE_SPLITTERS
            rm $MERGE_FULL $MERGE_DISCORDANTS $MERGE_SPLITTERS"
	fi
	$SAMBAMBA merge -t $THREADS $OUTPUT.bam $MERGE_FULL
	$SAMBAMBA merge -t $THREADS $OUTPUT.discordants.bam $MERGE_DISCORDANTS
	$SAMBAMBA merge -t $THREADS $OUTPUT.splitters.bam $MERGE_SPLITTERS
	rm $MERGE_FULL $MERGE_DISCORDANTS $MERGE_SPLITTERS

	# index the files
	if [[ $VERBOSE -eq 1 ]]
	then
	    echo -e "
	    $SAMBAMBA index $OUTPUT.bam
	    $SAMBAMBA index $OUTPUT.discordants.bam
	    $SAMBAMBA index $OUTPUT.splitters.bam"
	fi
	echo "
	$SAMBAMBA index $OUTPUT.bam
	$SAMBAMBA index $OUTPUT.discordants.bam
	$SAMBAMBA index $OUTPUT.splitters.bam
	" | $PARALLEL -j 3

    fi

    # clean up
    rm -r $TEMP_DIR

    echo "Done"

    # exit cleanly
    exit 0
}

# Show usage when there are no arguments.
if test -z "$1"
then
    usage
    exit 1
fi

while getopts "K:h" OPTION
do
    case $OPTION in
        h)
            usage
            exit 1
            ;;
	K)
	    ;;
        ?)
            usage
            exit
            ;;
    esac
done

# call the function
case "$1" in 
    'align')
	align "${@:2}"
	;;
    'var')
	var "${@:2}"
	;;
    'somatic')
	somatic "${@:2}"
	;;
    'sv')
	sv "${@:2}"
	;;
    'realign')
	realign "${@:2}"
	;;
    *)
	usage
	echo -e "Error: command \"$1\" not recognized\n"
	exit 1
esac

## END SCRIPT