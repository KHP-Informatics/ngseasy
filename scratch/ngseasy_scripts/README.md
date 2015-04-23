NGSeasy
===================
**A [Dockerized](https://www.docker.com/) and [Virtulaized](https://www.virtualbox.org/) ngs pipeline and tool-box.** 

**With NGSeasy you can now have full suite of NGS tools up and running on any high end workstation in an afternoon**

**Note:** NGSeasy is under **continous development** and the dev version evolves quickly. 

Please contact us for help/guidance on using the beat release. Our public server is down, so please to email us for access to the NGS resources and VM(s).  

- **NGSeasy-v1.0 Full Production release will be available Dec 2014**    

## Latest Version
*****************
[NGSeasy-v0.9.4](https://github.com/KHP-Informatics/ngs/tree/release-0.9.4) - NGSeasy: Latest Release  

- NB: so far full_no_gatk.sh is the most developed module. 

### Releases
[NGSeasy-v0.9.4](https://github.com/KHP-Informatics/ngs/tree/release-0.9.4) - NGSeasy: Latest Release  
[NGSeasy-v0.9.3](https://github.com/KHP-Informatics/ngs/tree/release-0.9.3) - NGSeasy: Minor Bug Fixes [full_no_gatk.sh]  
[NGSeasy-v0.9.2](https://github.com/KHP-Informatics/ngs/tree/release-0.9.2) - NGSeasy: Beta NovoIndex Fix  
[NGSeasy-v0.9.1](https://github.com/KHP-Informatics/ngs/tree/release-0.9.1) - NGSeasy: Beta added gatk cleaning steps  
[NGSeasy-v0.9](https://github.com/KHP-Informatics/ngs/tree/release-0.9)     - NGSeasy: Beta  

### Coming Soon
- Savant 
- SLOPE (CNV fo targetted NSG)  
- Cancer Pipelines
- Annotation Pipelines and Databases
- Visualisation Pipelines
- New Aligners:- GSNAP, mr- and mrs-Fast,gem
- Var Callers:- VarScan2
- SGE scripts and basic BASH scrips for running outside of Docker

### Authors
- Amos Folarin <amosfolarin@gmail.com>  
<a href="http://www.linkedin.com/pub/amos-folarin/34/b06/978">
<img src="http://www.linkedin.com/img/webpromo/btn_viewmy_160x33.png" width="160" height="33" alt="View Amos's profile on LinkedIn">
</a>

- Stephen J Newhouse <stephen.j.newhouse@gmail.com>  
<a href="http://uk.linkedin.com/pub/dr-stephen-newhouse/29/89a/11a">
<img src="http://www.linkedin.com/img/webpromo/btn_viewmy_160x33.png" width="160" height="33" alt="View Steve's profile on LinkedIn">
</a>

****************

We present **NGSeasy (Easy Analysis of Next Generation Sequencing)**, a flexible and easy-to-use NGS pipeline for automated alignment, quality control, variant calling and annotation. The pipeline allows users with minimal computational/bioinformatic skills to set up and run an NGS analysis on their own samples, in less than an afternoon, on any operating system (Windows, iOS or Linux) or infrastructure (workstation, cluster or cloud).

NGS pipelines typically utilize a large and varied range of software components and incur a substantial configuration burden during deployment which limits their portability to different computational environments. NGSeasy simplifies this by providing the pipeline components encapsulated in Docker™ containers and bundles in a wide choice of tools for each module. Each module of the pipeline represents one functional grouping of tools (e.g. sequence alignment, variant calling etc.).

Deploying the pipeline is as simple as pulling the container images from the public repository into any host running Docker. NGSeasy can be deployed on any medium to high-end workstation, high performance computer cluster and compute clouds (public/private cloud computing) - enabling instant access to elastic scalability without investment overheads for additional compute hardware and makes open and reproducible research straight forward for the greater scientific community.

- **NGSeasy updates every 6 months:**

**Lets us know if you want other tools added to NGSeasy**

***********

## Table of Contents
[NOTICE TO USERS](https://github.com/KHP-Informatics/ngs/blob/master/containerized/README.md#notice-to-users-of-the-container-image-or-vm)  

[Software requiring registration](https://github.com/KHP-Informatics/ngs/blob/master/containerized/README.md#software-composing-the-pipeline-requiring-registration)  

[Overview of Pipeline Components](https://github.com/KHP-Informatics/ngs/blob/master/containerized/README.md#overview-of-pipeline-components)  

[NGS Tools](https://github.com/KHP-Informatics/ngs/blob/master/containerized/README.md#the-tools-included-are-as-follows--)  

[Dockerised NGSEASY set up](https://github.com/KHP-Informatics/ngs/blob/master/containerized/README.md#dockerised-ngseasy)  

[Installing Docker](https://github.com/KHP-Informatics/ngs/blob/master/containerized/README.md#installing-docker)  

[Getting the Dockerised NGSEASY Pipeline](https://github.com/KHP-Informatics/ngs/blob/master/containerized/README.md#getting-the-dockerised-ngseasy-pipeline)  

[Running the Dockerised NGSEASY Pipeline](https://github.com/KHP-Informatics/ngs/blob/master/containerized/README.md#running-the-dockerised-ngseasy-pipeline)  

[NGSEASY-VM : An NGS Tool Box](https://github.com/KHP-Informatics/ngs/blob/master/containerized/README.md#ngseasy-vm--an-ngs-tool-box)  

[Installing Oracle VirtualBox](https://github.com/KHP-Informatics/ngs/blob/master/containerized/README.md#installing-oracle-virtualbox)  

[Getting the ngseasy-vm](https://github.com/KHP-Informatics/ngs/blob/master/containerized/README.md#getting-the-ngseasy-vm)  

[Installing the ngseasy-vm](https://github.com/KHP-Informatics/ngs/blob/master/containerized/README.md#installing-the-ngseasy-vm)  

******

## NOTICE TO USERS OF THE CONTAINER IMAGE OR VM

While the software used to build the image is composed of free software versions
some of the software has restrictions on use particularly for commercial 
purposes. Therefore if you wish to use this for commercial purposes, then you 
leagally have to approach the owners of the various components yourself!  

This pipeline uses a number of pieces of software which require registration. 
By using this you are agreeing to observe the Terms and Conditions of the 
relevant pieces of software that compose this pipeline.  

### Software composing the pipeline requiring registration

If you want to build the image from the Dockerfile then you need to get your 
own versions of (below) in the build directory:

   * novoalign http://www.novocraft.com/
   * Stampy http://www.well.ox.ac.uk/project-stampy
   * GATK https://www.broadinstitute.org/gatk/
   * ANNOVAR http://www.openbioinformatics.org/annovar/

[Back to The Begining](https://github.com/KHP-Informatics/ngs/blob/master/containerized/README.md#ngs-easy-v10)

******

Overview of Pipeline Components
================================
The basic pipeline contains all the basic tools needed for manipulation and 
quality control of raw fastq files (ILLUMINA focused), SAM/BAM manipulation,
alignment, cleaning (based on GATK best practises [http://www.broadinstitute.org/gatk/guide/best-practices]) and first pass
variant discovery. Separate containers are provided for indepth variant annotation,
structural variant calling, basic reporting and visualisations.  

![ngsEASY](figs/ngsEASY_component_visualisation.png "Dockerized NGS Pipeline")

# The Tools included are as follows :- 

### Fastq manipulation
- FASTQC
- SEQTK
- TRIMMOMATIC
- FASTX TOOLKIT

### Alignmnet
- BWA
- BOWTIE2
- STAMPY
- NOVOALIGN [FULL VERSION NOT AVAILABLE FOR PUBLIC USE]

### SAM/BAM Processing
- GATK
- PICARDTOOLS
- SAMTOOLS

### MISC
- BEDTOOLS
- VCFTOOLS
- BCFTOOLS

### VARIANT CALLING
- GATK
- SAMTOOLS/BCFTOOLS
- FREEBAYES
- PLATYPUS

### VARIANT ANNOTATION
- ANNOVAR
- SNPEFF
- VEP

### CNV/Structural Variant CALLING
- lumpy
- delly
- m-HMM
- cn.MOPS
- ExomeDepth  

**Software Versions**

- Trimmomatic-0.32
- bwa-0.7.10
- bowtie2-2.2.3
- novocraftV3.02.07.Linux3.0
- stampy-1.0.23
- samtools-0.1.19
- picard-tools-1.115
- GenomeAnalysisTK-3.2-2
- Platypus_0.7.4
- fastqc_v0.11.2  

[Back to The Begining](https://github.com/KHP-Informatics/ngs/blob/master/containerized/README.md#ngs-easy-v10)

******

Dockerised NGSeasy
==========================
![docker](figs/Docker_container_engine_logo.png "Docker")  

## Installing Docker

Follow the simple instructions in the links provided below  

- [Mac](https://docs.docker.com/installation/mac/)  
- [Windows](https://docs.docker.com/installation/windows/)
- [Ubuntu](https://docs.docker.com/installation/ubuntulinux/)

A full set of instructions for multiple operating systems are available on the [Docker website](https://docs.docker.com/installation/).

Getting the Dockerised NGSeasy Pipeline
-------------------------------------------

We have adapted the current best practices from the Genome Analysis Toolkit (GATK, http://www.broadinstitute.org/gatk/guide/best-practices)  for processing raw alignments in SAM/BAM format and variant calling. The current workflow, has been optimised for Illumina platforms, but can easily be adapted for other sequencing platforms, with minimal effort.  

As the containers themselves can be run as executables with pre-specified cpu and RAM resources, the orchestration of the pipeline can be placed under the control of conventional load balancers if this mode is required.  

*****************

### Available NGSeasy Docker images
Available to download at our [compbio Docker Hub](https://hub.docker.com/u/compbio)

### The Main Power Tool!
**ngseasy-alignment-public:v1.2**
Docker Hub: [compbio/ngseasy-alignment-public:v1.2](https://registry.hub.docker.com/u/compbio/ngseasy-alignment-public/)

This is a large image (4.989 GB) containing all the tools needed to go from raw ``.fastq`` files to aligned ``.BAM`` to SNP and small INDEL variant calls ``.vcf`` .

```bash
sudo docker pull compbio/ngseasy-alignment-public:v1.2
```

Getting All NGSeasy images, Reasources and Scripts
------------------------------

All Images can be pulled down from [Docker Hub](https://hub.docker.com/u/compbio/) using the script [get_NGSeasy.sh](https://github.com/KHP-Informatics/ngs/blob/master/containerized/get_NGSeasy.sh)


**NGSeasy Reasources**

This is a 25GB ``ngseasy_resources.tar.gz`` file containing :-  

- ``reference_genomes_b37.tgz`` b37 reference genomes indexed for use with all provided aligners (BWA, Bowtie2, Stampy, Novoalign) and annotation bed files for use with pipeline scripts
- ``gatk_resources.tar.gz`` gatk resources bundle
- ``fastq_example.tgz`` Example 75bp PE Illumina Whole Exome Sequence fastq data for **NA12878**
- Annotation Databases Coming in the next update 
 
**SFTP Login Details**

```
location: upload.brc.iop.kcl.ac.uk  
Port: 51515  
user: ngseasy  
password: ngseasy  
```

Example:- 
```bash
$ sftp  ngseasy@upload.brc.iop.kcl.ac.uk
ngseasy@upload.brc.iop.kcl.ac.uk's password: 
Connected to upload.brc.iop.kcl.ac.uk.
sftp> ls
NGSeasy  
sftp> cd NGSeasy
sftp> ls
ngseasy-vm-v1.0.vdi          ngseasy_resources.tar.gz     
sftp> get -r *
```
I would recommend using a separate program like [FileZilla](https://filezilla-project.org/), which will make it much easier for you to set up and manage your file transfers

**NGSeasy pipeline scripts**


Clone latest NGSeasy scripts from out GitHub repository
```sh
git clone https://github.com/KHP-Informatics/ngs.git
```

[Back to The Begining](https://github.com/KHP-Informatics/ngs/blob/master/containerized/README.md#ngs-easy-v10)

****************

NGSeasy Project Set up
--------------------------

## Step 0. Set up project directories

Fastq files must have suffix and be gzipped: **_1.fq.gz** or **_2.fq.gz**  
furture version will allow any format  

```{bash}
#--------------------------------#
# make top level dirs 
#--------------------------------#
cd media
mkdir ngs_projects
mkdir ngs_projects/fastq_raw
mkdir ngs_projects/config_files
```
****

## Step 1. Set up project configuration file

In Excel make config file and save as [TAB] Delimited file with ``.tsv`` extenstion.  
See Example provided and [GoogleDoc](https://docs.google.com/spreadsheets/d/1kp1Nyw0x3zXqO2Wm2Z25ErJ0Z-Uoab8tjRPq9h4sonk/edit?usp=sharing). Remove the header from this file before running the pipeline. This sets up Information related to: Project Name, Sample Name, Library Type, Pipeline to call, NCPU.

The [config.file.tsv] should contain the following 15 columns for each sample to be run through a pipeline:- 

|Variable|type|Description|
|--------|--------|--------|
POJECT_ID|string|Project ID|
SAMPLE_ID|string|Sample ID|
FASTQ1|string|Raw fastq file name read 1|
FASTQ2|string|Raw fastq file name read 1|
PROJECT_DIR|string|Project Directory|
DNA_PREP_LIBRARY_ID|string|DNA Libray Prep ID|
NGS_PLATFORM|string|Platform Name|
NGS_TYPE|string|Experiment type|
BED_ANNO|string|Annotation Bed File
PIPELINE|string|NGSeasy Pipeline Script|
ALIGNER|string|Aligner|
VARCALLER|string|Variant Caller|
GTMODEGATK|string|GATK Variant Caller Mode|
CLEANUP|string|Clean Up Files (TRUE/FALSE)|
NCPU|number|Number of cores to call|

****

## Step 2. Initiate the Project
The user needs to make the relevent directory structure on their local machine. Running this script ensures that all relevant directories are set up, ans also enforces a clean structure to the NGS project.  

On our sysetm we typically set up a top-level driectory called `ngs_projects` within which we store output from all our individual NGS projects. Within this we make a `raw_fastq` folder, where we temporarily store all the raw fastq files for each project. This folder acts as an initial stagging area for the raw fastq files. During the project set up, we copy/move project/sample related fastq files to their own specific directories.

Running `ngseasy_initiate_project` with the relevent configuration file, will set up the following directory structure for every project and sample within a project:-  

### NGS Project Directory 
```
.
ngs_projects  
|  
|__raw_fastq  
|__config_files  
|__reference_genomes_b37  
|__gatk_resources  
|__ngseasy
|
|__ project_id  
	|  
	|__run_logs  
	|__config_files  
	|__project_vcfs  
	|__project_bams  
	|__project_reports  
	|
	|__sample_id  
		|  
		|__fastq  
		|__tmp  
		|__alignments  
		|__vcf  
		|__reports  
		|__config_files  

```
### Running **ngseasy_initiate_project**

```{bash}
ngseasy_initiate_project -c config.file.tsv -d /media/ngs_projects
```
****

## Step 3. Copy Project Fastq files to relevent Project/Sample Directories

```{bash}
ngseasy_initiate_fastq -c config.file.tsv -d /media/ngs_projects
```

****

## Step 4. Start the NGSeasy Volume Contaier

In the Docker container the project directory is mounted in `/home/pipeman/ngs_projects`

```{bash}
ngseasy_volumes_container -d /media/ngs_projects
```

## Summary

```{bash}
#make top level dirs 

mkdir ngs_projects
mkdir ngs_projects/fastq_raw
mkdir ngs_projects/config_files

#copy/download raw fastq file to [ngs_projects/fastq_raw]

#set up project specific configuration file [config.file.tsv]

# Set up NGSeasy
ngseasy_initiate_project -c config.file.tsv -d /media/ngs_projects
ngseasy_initiate_fastq -c config.file.tsv -d /media/ngs_projects
ngseasy_volumes_container -d /media/ngs_projects
```

** You are now ready to run **

****

## Running a full pipeline : from raw fastq to vcf calls

See https://github.com/KHP-Informatics/ngs/tree/dev2/bin for dev functions (Still workin on these). Each of these will call a separate container and run a part of the NGS pipeline. Each step is usually dependent on the previous step(s) - in that they require certain data/input/output in the correct format and with the correct nameing conventions enforced by our pipeline to exist, before executing.

A full pipeline is set out below :-  

```{bash}

#--------------------------------#
# make top level dirs 
#--------------------------------#
cd media
mkdir ngs_projects
mkdir ngs_projects/fastq_raw # fastq staging area
mkdir ngs_projects/config_files # config files
mkdir ngs_projects/humandb # for annovar databses

#get NGSeasy resources
# sftp From ........copy data to and extract
cd ngs_projects
sftp....

tar xvf gatk_resources.tar.gz; gunzip *
tar xvf reference_genomes_b37.tgz; gunzip *


#--------------------------------#
# get and PATH nsgeasy scripts
#--------------------------------#

cd ngs_projects/nsgeasy
git clone https://github.com/KHP-Informatics/ngs.git
git checkout dev2

# eg :- 

export PATH=$PATH:/media/ngs_projects/nsgeasy/ngs/bin

# or add to global .bashrc

echo "export PATH=$PATH:/media/ngs_projects/nsgeasy/ngs/bin" ~/.bashrc
source ~/.bashrc


#to do [get_annovar_humandb]

#get images
bash get_containers.sh

#------------------------------------------------#
# to be run outside of docker and before ngseasy #
#------------------------------------------------#

ngseasy_initiate_project -c config.file.tsv -d /media/ngs_projects

ngseasy_initiate_fastq -c config.file.tsv -d /media/ngs_projects

ngseasy_get_annovar_db -d /media/ngs_projects/humandb

#--------------------#
# NGSEASY Dockerised #
#--------------------#

# A pipeline is called using :-

    ngseasy -c config.file.tsv -d /media/nsg_projects

# in the config file we as to call the pipeline [full]
# here [ngs_full_gatk] is a wrapper/fucntion for calling the pipeline

Each function is a bash wrapper for an image/container(s)



ngs_full_gatk() { 

#step through a full pipeline

#1 FastQC raw data
ngseasy_fastqc -c ${config_tsv} -d ${project_directory} -t 1;

#2 Trimm Fastq files
ngseasy_fastq_trimm -c ${config_tsv} -d ${project_directory};

#3 FastQC trimmed file
ngseasy_fastqqc -c ${config_tsv} -d ${project_directory} -t 0;

#4 Run Alignment
ngseasy_alignment -c ${config_tsv} -d ${project_directory};

#5 Add read groups ie sample info to BAM file
nsgeasy_add_readgroups -c ${config_tsv} -d ${project_directory};

#5 Mark Duplicates
ngseasy_mark_dupes -c ${config_tsv} -d ${project_directory};

#6 Indel Realignment
nsgeasy_indel_realignment -c ${config_tsv} -d ${project_directory};

#7 Base Recalibration
nsgeasy_base_recalibration -c ${config_tsv} -d ${project_directory};

#8 Alignment QC reports
ngseasy_alignment_qc -c ${config_tsv} -d ${project_directory};

#9 Call Vairants : SNPS and small Indels
nsgeasy_var_callers -c ${config_tsv} -d ${project_directory};

#10 CNV Calling
nsgeasy_cnv_callers -c ${config_tsv} -d ${project_directory};

#11 Variant Annotation: Annovar
nsgeasy_variant_annotation -c ${config_tsv} -d ${project_directory};

#TODO
#variant filter
#nsgeasy_variant_combine -c ${config_tsv} -d ${project_directory}
#nsgeasy_variant_report -c ${config_tsv} -d ${project_directory}


}

recommend full : trimmed aln gatk filtered and ensemble calls (multi SNP/INDELS/SV callers) base recalibration
if not novoalign then stampy (bwa with stampy)


# Pipelines :  


ngs_full_gatk  
ngs_full_no_gatk  
ngs_fq2bam # 





```
****

To add :
- getting the pipeline and setting up resources data
- not all steps need config?
- pipeline option need to be set how? list of steps, specified full, full_no_gatk, var_call_only, cnv_call_only, qc_reports??

****

Output suffixes 
===================

### Alignment Output
*.raw.sam  (WEX ~ 8GB)
*.raw.bam  
*.raw.bai  
*.sort.bam  (WEX ~ 3GB)
*.sort.bai  

****
### Addreadgroup
*.addrg.bam  
*.addrg.bai  
*.addrg.bam.bai  

****
### Dupemark
*.dupemk.bam  
*.dupemk.bai
*.dupemk.bam.bai  

***
### Indel realign
*.realn.bam  
*.realn.bai
*.realn.bam.bai  


***
### Base recal
*.recal.bam  (WEX ~ 4.4G)  
*.recal.bai  
*.recal.bam.bai  
*.realn.bam.BaseRecalibrator.table  
*.recal.bam.BaseRecalibrator.table  
*.recal.bam.BaseRecalibrator.BQSR.csv  

***

## Thresholds for Variant calling etc

For Freebayes and Platypus tools:-  

- We set min coverage to 10  
- Min mappinng quality to 20  
- Min base quality to 20

For GATK HaplotypeCaller (and UnifiedGenotyper)

```-stand_call_conf 30 -stand_emit_conf 10 -dcov 250 -minPruning 10```

Note: ```minPruning 10``` was added as many runs of HaplotypeCaller failed when using non-bwa aligend and GATK best practices cleaned BAMs. This fix sorted all problems out, and you really dont want dodgy variant calls...do you? Same goes for thresholds hard coded for use with Freebayes and Platypus.  
These setting all work well in our hands. Feel  free to edit the scripts to suit your needs.


****
blah blah blah

****

## Gottchas

**bin/bash -c**

- need to add ```/bin/bash -c ${COMMAND}``` when software require ```>``` redirect to some output

example below for bwa:-  

```
  sudo docker run \
  -P \
  --name sam2bam_${SAMPLE_ID} \
  --volumes-from volumes_container \
  -t compbio/ngseasy-samtools:v0.9 /bin/bash -c \
  "/usr/local/pipeline/samtools/samtools view -bhS ${SOUTDocker}/alignments/${BAM_PREFIX}.raw.bwa.sam > ${SOUTDocker}/alignments/${BAM_PREFIX}.raw.bwa.bam"
  ```

runnig this without ```/bin/bash -c``` breaks. The ```>``` is called outside of the container

### The Annoying thing about GATK!
This will break your runs if multiple calls try and access the file when the first call deletes it!  
```
WARN  11:05:27,577 RMDTrackBuilder - Index file /home/pipeman/gatk_resources/Mills_and_1000G_gold_standard.indels.b37.vcf.idx is out of date (index older than input file), deleting and updating the index file 
INFO  11:05:31,699 RMDTrackBuilder - Writing Tribble index to disk for file /home/pipeman/gatk_resources/Mills_and_1000G_gold_standard.indels.b37.vcf.idx 
```


## CNV tools to think about
EXCAVATOR: detecting copy number variants from whole-exome sequencing data @ http://genomebiology.com/2013/14/10/R120


>We developed a novel software tool, EXCAVATOR, for the detection of copy number variants (CNVs) from whole-exome sequencing data. EXCAVATOR combines a three-step normalization procedure with a novel heterogeneous hidden Markov model algorithm and a calling method that classifies genomic regions into five copy number states. We validate EXCAVATOR on three datasets and compare the results with three other methods. These analyses show that EXCAVATOR outperforms the other methods and is therefore a valuable tool for the investigation of CNVs in largescale projects, as well as in clinical research and diagnostics. EXCAVATOR is freely available at http://sourceforge.net/projects/excavatortool/ webcite.

## AnalyzeCovariates post recal plots ~ needs R rehsape and gplots

AnalyzeCovariates keeps failing - due to R. GATK dont document R requirements for these scripts to work!  

```
library("ggplot2")
library(gplots)
library("reshape")
library("grid")
library("tools") #For compactPDF in R 2.13+
library(gsalib)
```

URL: https://github.com/broadgsa/gatk-protected/blob/323f22f852c90a6cf53ece5e72c165b1166ad8c7/public/gatk-tools-public/src/main/resources/org/broadinstitute/gatk/utils/recalibration/BQSR.R

```
install.packages(c("ggplot2","gplots","reshape","grid","tools","gsalib"),dependencies=TRUE)

Warning message:
packages 'grid', 'tools' are not available (for R version 3.1.1) 
```
Drop plots of recal data

****

## bcbio.variation
https://github.com/chapmanb/bcbio.variation

# bcbio.variation

A toolkit to analyze genome variation data, built on top of the
[Genome Analysis Toolkit (GATK)][1] with Clojure. It supports scoring for the
[Archon Genomics X PRIZE competition][5] and is also a general framework for
variant file comparison. It enables validation of variants and exploration of
algorithm differences between calling methods by automating the process involved
with comparing two sets of variants. For users, this integrates with the
[bcbio-nextgen][8] framework to automate variant calling and validation. For
developers, bcbio.variation provides command line tools and an API to clean and
normalize variant data in [VCF format][2] avoiding comparison artifacts
associated with different variant representations.

* [Description of the comparison framework and variant calling algorithm comparisons][7]
* [Presentation from Bioinformatics Open Source Conference 2012][p1]
* [Presentation overview of the project][4]
* [Howto description of interfacing with GATK][6]
* [Code documentation][3]

[![Build Status](https://secure.travis-ci.org/chapmanb/bcbio.variation.png)](http://travis-ci.org/chapmanb/bcbio.variation)

[1]: http://www.broadinstitute.org/gsa/wiki/index.php/The_Genome_Analysis_Toolkit
[2]: http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40
[3]: http://chapmanb.github.com/bcbio.variation
[4]: http://chapmanb.github.com/bcbio.variation/presentations/gatk_clojure.pdf
[5]: http://genomics.xprize.org/
[6]: http://bcbio.wordpress.com/2012/03/04/extending-the-gatk-for-custom-variant-comparisons-using-clojure/
[p1]: http://chapmanb.github.com/bcbio.variation/presentations/variation_bosc_2012/variation_chapman.pdf
[7]: http://bcbio.wordpress.com/2013/05/06/framework-for-evaluating-variant-detection-methods-comparison-of-aligners-and-callers/
[8]: https://github.com/chapmanb/bcbio-nextgen

******

## Java issues in docker 
All java set as 1.6. Updated and commit in docker images interactively

```
sudo update-alternatives --config java
```

****


#biobambam
https://github.com/gt1/biobambam

biobambam
======

This package contains some tools for processing BAM files including

 - bamcollate2: reads BAM and writes BAM reordered such that alignment
   or collated by query name
 - bammarkduplicates: reads BAM and writes BAM with duplicate alignments
   marked using the BAM flags field
 - bammaskflags: reads BAM and writes BAM while masking (removing)
   bits from the flags column
 - bamrecompress: reads BAM and writes BAM with a defined compression
   setting. This tool is capable of multi-threading.
 - bamsort: reads BAM and writes BAM resorted by coordinates or query
   name
 - bamtofastq: reads BAM and writes FastQ; output can be collated or
   uncollated by query name

A short list of options is available for each program by calling it
with the -h parameter, e.g.

	bamsort -h

Source
------

The biobambam source code is hosted on github:

	git@github.com:gt1/biobambam.git

Compilation of biobambam
------------------------

biobambam needs libmaus [https://github.com/gt1/libmaus] . When libmaus
is installed in ${LIBMAUSPREFIX} then biobambam can be compiled and
installed in ${HOME}/biobambam using

	- autoreconf -i -f
	- ./configure --with-libmaus=${LIBMAUSPREFIX} \
		--prefix=${HOME}/biobambam
	- make install

	***