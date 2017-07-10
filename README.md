![NGSeasy_logo](https://github.com/KHP-Informatics/ngseasy/blob/sjn_dev/figs/NGSeasy_logo_0.0.1.png)

		
		I am looking for collaborators (dev help). Please get in contact

NGSeasy: A [Dockerized](https://www.docker.com/) NGS pipeline and tool-box
=====================================================================================
[![Join the chat at https://gitter.im/KHP-Informatics/ngseasy](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/KHP-Informatics/ngseasy?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge) [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.31444.svg)](http://dx.doi.org/10.5281/zenodo.31444)

## NGSeasy

**With NGSeasy you can now have full suite of NGS tools up and running on any high end workstation in an afternoon**  

We present NGSeasy (Easy Analysis of Next Generation Sequencing), a flexible and easy-to-use NGS pipeline for automated alignment, quality control, variant calling and annotation. The pipeline allows users with minimal computational/bioinformatic skills to set up and run an NGS analysis on their own samples, in less than an afternoon, on any operating system (Windows, iOS or Linux) or infrastructure (workstation, cluster or cloud).

**Authors:** Stephen J Newhouse and Amos Folarin  
**Release Version:** **1.0-r001**  
**Release:** dirty_tango  
**Citation:** [Folarin AA, Dobson RJ and Newhouse SJ. NGSeasy: a next generation sequencing pipeline in Docker containers [version 1; referees: 3 approved with reservations] F1000Research 2015, 4(ISCB Comm J):997 (doi: 10.12688/f1000research.7104.1).](http://f1000research.com/articles/4-997/v1)   

## Acknowledgements
This work is funded by the National Institute for Health Research (NIHR) Biomedical Research Centre at South London and Maudsley NHS Foundation Trust and King’s College London.

- **Lets us know if you want other tools added to NGSeasy**

[**NGSeasy is completely open source and we encourage interested folks to jump in and get involved in the dev with us.**](https://github.com/KHP-Informatics/ngseasy.git)

## Contributing to NGSEASY
- Fork it!
- Create your **feature branch**: `git checkout -b my-new-feature`
- Commit your changes: `git commit -am 'Add some feature'`
- Push to the branch: `git push origin my-new-feature`
- Submit a pull request!

## NGSeasy: Genome Comparison & Analytic Testing ([GCAT](http://www.bioplanet.com/gcat)) Reports
Here we provide a quick look at basic NGSeasy performance (more results coming soon).

|GCAT Report|Test Data |Pipeline|
|----|----|----|
|[**NGSEASY-NTRIM-BWA-FREEBAYES-D**](http://www.bioplanet.com/gcat/reports/6167-seeirhwtfp/variant-calls/illumina-100bp-pe-exome-150x/ngseasy-ntrim-bwa-freebayes-d/compare-570-270-181/group-read-depth)|[illumina-100bp-pe-exome-150x](https://s3-eu-west-1.amazonaws.com/ngseasy.data/fastq_test_data/illumina.100bp.pe.wex.150x_1.fastq.gz)|[fastq > bwa > freebayes](https://github.com/KHP-Informatics/ngseasy/blob/master/test/ngseasy_test.config.freebayes.tsv)|
|[**NGSEASY-NTRIM-BWA-PLATYPUS-D**](http://www.bioplanet.com/gcat/reports/6264-uzdswfzmkv/variant-calls/illumina-100bp-pe-exome-150x/ngseasy-ntrim-bwa-platypus-d/compare-571-180/group-read-depth)|[illumina-100bp-pe-exome-150x](https://s3-eu-west-1.amazonaws.com/ngseasy.data/fastq_test_data/illumina.100bp.pe.wex.150x_2.fastq.gz)|[fastq > bwa > platypus](https://github.com/KHP-Informatics/ngseasy/blob/master/test/ngseasy_test.config.platypus.tsv)|

An example of the run commands:

```bash
ngseasy -c ngseasy_test.config.freebayes.tsv -d /media/Data/ngs_projects
ngseasy -c ngseasy_test.config.platypus.tsv  -d /media/Data/ngs_projects
```

## Author Contact Details

Please contact us for help/guidance on using the beta release.

|Author|email|Twitter| LinkedIn|
|----|----|----|----|
|Dr Stephen J Newhouse| <stephen.j.newhouse@gmail.com>| [@s_j_newhouse](https://twitter.com/s_j_newhouse?lang=en)| View Steve's profile on [LinkedIn](http://uk.linkedin.com/pub/dr-stephen-newhouse/29/89a/11a)|
|Dr Amos Folarin| <amosfolarin@gmail.com>| [@amosfolarin](https://twitter.com/amosfolarin?lang=en)|View Amos's profile on [LinkedIn](http://www.linkedin.com/pub/amos-folarin/34/b06/978)|



## Issues, Questions and Queries

**Please Direct all queries to [https://github.com/KHP-Informatics/ngseasy/issues]**

When sending bug reports etc please provide:-

- Date of Download
- OS and version
- Basic Machine Specs (CPU, RAM)
- Docker Version  
- Network Speed ([Testing Internet Connection Speed](https://github.com/KHP-Informatics/ngseasy/tree/master#recommended-network-speed))
- The Code you ran eg:- `ngseasy -c my.config.tsv -d /My/Dir `
- your config file  
- The exact error as printed to screen

**WARNING!**
NGSeasy is not [numpty](http://www.urbandictionary.com/define.php?term=Numpty) or bad data proof!  

**Please read the docs**, stay calm, take your time and think about what you are doing...and if [www.google.com] doesnt help, then please direct all queries to [https://github.com/KHP-Informatics/ngseasy/issues].

## Install Docker

Full instructions at https://docs.docker.com/.

Some fixes to make life easy...allows you to run ``docker`` without ``sudo``.

This may differ for your OS, and mostly applies to flavours of ``Linux``. Check with your sys admin or just Google https://www.google.com.

MAC/Windows users using http://boot2docker.io/ should be fine. Read the docs or just Google https://www.google.com.

**Create a docker group**

```bash
sudo addgroup docker
```

**Add user to docker group**  

Here user is ``ec2-user``

```bash
sudo usermod -aG docker ec2-user
```

Log out and log back in.

This ensures your user is running with the correct permissions.

Verify your work by running ``docker`` without ``sudo``.

```bash
docker run hello-world
```

..this is what you should get...

```
Unable to find image 'hello-world:latest' locally
Pulling repository hello-world
91c95931e552: Download complete
a8219747be10: Download complete
Status: Downloaded newer image for hello-world:latest
Hello from Docker.
This message shows that your installation appears to be working correctly.

To generate this message, Docker took the following steps:
 1. The Docker client contacted the Docker daemon.
 2. The Docker daemon pulled the "hello-world" image from the Docker Hub.
    (Assuming it was not already locally available.)
 3. The Docker daemon created a new container from that image which runs the
    executable that produces the output you are currently reading.
 4. The Docker daemon streamed that output to the Docker client, which sent it
    to your terminal.

To try something more ambitious, you can run an Ubuntu container with:
 $ docker run -it ubuntu bash

For more examples and ideas, visit:
 http://docs.docker.com/userguide/
```

## Docker Security...

This post reviews the various security implications of using Docker to run applications within containers, and how to address them:
[How Secure are Containers?](https://blog.docker.com/2013/08/containers-docker-how-secure-are-they/#more-697)

> Docker containers are, by default, quite secure; especially if you take care of running your processes inside the containers as non-privileged users (i.e. non root).

## Get NGSeasy

```bash

#############################################
## Get NGSeasy                             ##
#############################################

cd /home/${USER}

git clone https://github.com/KHP-Informatics/ngseasy.git
```

## Install NGSeasy

- Default install directory is `/home/${USER}`
- in this example user home is `/home/ec2-user`
- `make INSTALLDIR="/home/ec2-user" all`
   - sets up top level directory structure
   - gets all docker images
   - gets indexed hg19 and b37 genomes
   - gets GATK recources for hg19 and b37 genomes
   - gets whole genome and exome test data
- Always set your **`INSTALLDIR`** : If you run `sudo make all` the install path will be `/home/root`. Please dont do this!
- `sudo make install` installs scripts to `/usr/local/bin/`  

```bash
#############################################
## install NGSeasy                         ##
#############################################

cd ngseasy

## 1.
make INSTALLDIR="/home/ec2-user" all

## 2.
sudo make install
```

Installation can take a while, 1-2 hours, so go get a coffee../just chill...if your network is bad...then who knows how long...still..just chill...or go get fast internet!

## NGSeasy Security
**All NGSeasy applications are run as the non-root user `pipeman` within each container**

## Recommended Network Speed

**> 500 Mbit/s** : anything less will add a lot of time to set up (days - weeks).

**Testing Internet Connection Speed**

source : http://askubuntu.com/questions/104755/how-to-check-internet-speed-via-terminal

```bash
wget -O speedtest-cli https://raw.github.com/sivel/speedtest-cli/master/speedtest_cli.py
chmod +x speedtest-cli
./speedtest-cli
```

```
Retrieving speedtest.net configuration...
Retrieving speedtest.net server list...
Testing from Comcast Cable (x.x.x.x)...
Selecting best server based on ping...
Hosted by FiberCloud, Inc (Seattle, WA) [12.03 km]: 44.028 ms
Testing download speed........................................
Download: 32.29 Mbit/s
Testing upload speed..................................................
Upload: 5.18 Mbit/s
```

## Install time on Amazon EC2

Connection Speed: ~ 800 Mbit/s

```
real    94m54.237s
user    12m26.960s
sys     28m46.648s
```

Note: We have only tested NGSeasy installation on Amazon EC2, Openstack and UK University Networks. These are all fairly fast networks with speeds exceeding 800 Mbit/s on average.

## Running NGSeasy for the first time on the test data

**Important!** NGSeasy is controlled from a single `config` file. See [ngseasy_test.config.tsv](https://github.com/KHP-Informatics/ngseasy/blob/master/test/ngseasy_test.config.tsv) for a basic template. It is important that the user sets this up properly before running NGSeasy.

```bash
#############################################
## 0. Move to config file dir

cd /home/ec2-user/ngs_projects/config_files/

#############################################
## 1. Run basic test

ngseasy -c ngseasy_test.config.tsv -d /home/ec2-user/ngs_projects
```

### What should happen...

This runs the following basic pipeline on Whole Exome PE 30x Illumina data, aligning to b37 (in theory...give it a try).  

- FastQC > Trimmomatic > BWA > Platypus

### Some notes and pointers

- Edit **NCPU** in  **[ngseasy_test.config.tsv]** to suit your system  
- Edit **PROJECT\_DIR** in  **[ngseasy_test.config.tsv]** to suit your install path
- We expect the user to palce all raw fastq files in ``raw_fastq``. NGSeasy uses this
as a stagging area for new project and sample data.  
- right now, always run ``ngseasy`` from the location/directory that contains the config.file  
- each component of ngseasy can be run as a standalone script  

****************


NGSeasy (Easy Analysis of Next Generation Sequencing)
=======================================================
We present **NGSeasy (Easy Analysis of Next Generation Sequencing)**, a flexible and easy-to-use NGS pipeline for automated alignment, quality control, variant calling and annotation. The pipeline allows users with minimal computational/bioinformatic skills to set up and run an NGS analysis on their own samples, in less than an afternoon, on any operating system (Windows, iOS or Linux) or infrastructure (workstation, cluster or cloud).

NGS pipelines typically utilize a large and varied range of software components and incur a substantial configuration burden during deployment which limits their portability to different computational environments. NGSeasy simplifies this by providing the pipeline components encapsulated in Docker™ containers and bundles in a wide choice of tools for each module. Each module of the pipeline represents one functional grouping of tools (e.g. sequence alignment, variant calling etc.).

Deploying the pipeline is as simple as pulling the container images from the public repository into any host running Docker. NGSeasy can be deployed on any medium to high-end workstation, high performance computer cluster and compute clouds (public/private cloud computing) - enabling instant access to elastic scalability without investment overheads for additional compute hardware and makes open and reproducible research straight forward for the greater scientific community.

### Advantages ###
- Easy to use for non-informaticians.  
- All run from a single config file that can be made in Excel.  
- User can select from mutiple aligners, variant callers and variant annotators
- No scary python, .yaml or .json files...just one simple Excel workbook saved as a textfile.  
- Just follow our simple set of instructions and NGS away!  
- Choice of aligners and variant callers and anntators  
- Allows reproducible research  
- Version controlled for auditing  
- Customisable  
- Easy to add new tools  
- If it's broke...we will fix it..
- Enforced naming convention and directory structures  
- Allows users to run "Bake Offs" between tools with ease  

We have adapted the current best practices from the Genome Analysis Toolkit (GATK, http://www.broadinstitute.org/gatk/guide/best-practices)  for processing raw alignments in SAM/BAM format and variant calling. The current workflow, has been optimised for Illumina platforms, but can easily be adapted for other sequencing platforms, with minimal effort.  

As the containers themselves can be run as executables with pre-specified cpu and RAM resources, the orchestration of the pipeline can be placed under the control of conventional load balancers if this mode is required.  

## Genomes

| Genome Build | Status |
|--------|------|
|hs37d5| coming soon |
|b37|available|
|hg19|available|
|hs38DH|coming soon |

****************

Overview of the NGSeasy Pipeline Components
=============================================
The basic pipeline contains all the basic tools needed for manipulation and
quality control of raw fastq files (ILLUMINA focused), SAM/BAM manipulation,
alignment, cleaning (based on GATK best practises [http://www.broadinstitute.org/gatk/guide/best-practices]) and first pass
variant discovery. Separate containers are provided for indepth variant annotation,
structural variant calling, basic reporting and visualisations.  

![ngsEASY](https://github.com/KHP-Informatics/ngs/blob/dev2/figs/ngsEASY_component_visualisation.png)

## A Special note on the NGSeasy base image.

We include the following - what we think of as - **_NGS Powertools_** in the **[compbio/ngseasy-base]()** image.
These are all tools that allow the user to slice and dice BED/SAM/BAM/VCF files in multiple ways.

 1.  samtools  
 2.  bcftools  
 3.  vcftools  
 4.  vcflib  
 5.  bamUtil  
 6.  bedtools2  
 7.  ogap  
 8.  samblaster  
 9.  sambamba  
 10. bamleftalign  
 11. seqtk  
 12. parallel  

This image is used as the base of all our compbio/ngseasy-* tools.

**Why not a separate containers per application?** The more docker-esque approach, would be to have separate containers for each NGS tool. However, this belies the fact that many of these tools interact in a deep way. Therefore, we built  these into a single development environment for ngseasy, to allow pipes and streamlined system calls for manipulating the output of NGS pipelines (BED/SAM/BAM/VCF files).

****************

The Full NGSeasy pipeline
=============================

The NGSeasy pipelines implement the following :-

- **Quality control of raw fastq** files using **[FASTQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)**  

- **Read trimming** using **[TRIMMOMATIC](http://www.usadellab.org/cms/?page=trimmomatic)**.

- **Alignment** using one of
    - **[BWA](http://bio-bwa.sourceforge.net/)**  
    - **[STAMPY](http://www.well.ox.ac.uk/project-stampy)**
    - **[NOVOALIGN](http://www.novocraft.com)**  
    - **[BOWTIE2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)**  
    - **[SNAP](http://snap.cs.berkeley.edu/)**  

- **SAM/BAM sorting and indexing** with **[SAMBAMBA](https://github.com/lomereiter/sambamba)**.  

- **Read Group information added** using **[PICARDTOOLS](http://broadinstitute.github.io/picard/):[AddOrReplaceReadGroups](http://broadinstitute.github.io/picard/command-line-overview.html#AddOrReplaceReadGroups)**

- **Duplicate marking** with **[SAMBLASTER](https://github.com/GregoryFaust/samblaster)**.  

>For academic users and/or commercial/clinical groups whom have paid for GATK licensing, the next steps are to perform

- **Indel indel realignment and base quality score recalibration** using **[GATK](https://www.broadinstitute.org/gatk/)** built in tools :
    - **[GATK:RealignerTargetCreator](https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_indels_RealignerTargetCreator.php)**
    - **[GATK:IndelRealigner](https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_indels_IndelRealigner.php)**
    - **[GATK:BaseRecalibrator](https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_bqsr_BaseRecalibrator.php)**

> For the non-GATK version

- **Base quality score recalibration** using [BamUtil](http://genome.sph.umich.edu/wiki/BamUtil)  
    - **[BamUtil:recab](http://genome.sph.umich.edu/wiki/BamUtil:_recab)**

- **Post alignment quality control and reporting** is performed usng a number of tools and custom scripts:
    - **[SAMTOOLS:flagstats](https://github.com/samtools/samtools)**
    - **[BEDTOOLS:genomecov](https://github.com/arq5x/bedtools2)**
    - **[BEDTOOLS:bamtobed](https://github.com/arq5x/bedtools2)**
    - **[PICARDTOOLS:CollectMultipleMetrics](http://broadinstitute.github.io/picard/command-line-overview.html#CollectMultipleMetrics)**
    - **[PICARDTOOLS:CollectAlignmentSummaryMetrics](http://broadinstitute.github.io/picard/command-line-overview.html#CollectAlignmentSummaryMetrics)**
    - **[PICARDTOOLS:CollectWgsMetrics](http://broadinstitute.github.io/picard/command-line-overview.html#CollectWgsMetrics)**
    - **[PICARDTOOLS:CollectTargetedPcrMetrics](http://broadinstitute.github.io/picard/command-line-overview.html#CollectTargetedPcrMetrics)** (coming soon)

- **SNP and small INDEL** calling using one of the following or a combibation of these tools, if the `ensemble` method is called using **[bcbio.variation variant-ensemble](https://github.com/chapmanb/bcbio.variation)**
    - **[FREEBAYES](https://github.com/ekg/freebayes)**
    - **[PLATYPUS](http://www.well.ox.ac.uk/platypus)**
    - **[GATK:UnifiedGenotyper](https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_genotyper_UnifiedGenotyper.php)**
    - **[GATK:HaplotypeCaller](https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php)**

- **Structural Variant (CNV)** calling using one of the following or or a combibation of if the `ensemble` methods are called:-
    - **[DELLY](https://github.com/tobiasrausch/delly) : in Dev**
    - **[LUMPY](https://github.com/arq5x/lumpy-sv/): in Dev**
    - **[cn.MOPS](http://www.bioinf.jku.at/software/cnmops/): in Dev**
    - **[m-HMM](https://www.stt.msu.edu/users/hengwang/mHMM.html): in Dev**
    - **[ExomeDepth](http://cran.r-project.org/web/packages/ExomeDepth/index.html): in Dev**
    - **[SLOPE](http://www-genepi.med.utah.edu/suppl/SLOPE/index.html): in Dev and not tested**  
    - **[cnvnator](http://sv.gersteinlab.org/): in Dev**  

- **Variant annotation** using using one of the following or a combibation of if the `ensemble` methods are called.
    - **[SnpEff](http://snpeff.sourceforge.net/): in Dev**
    - **[ANNOVAR](http://www.openbioinformatics.org/annovar/): in Dev**
    - **[VEP](http://www.ensembl.org/info/docs/tools/vep/index.html): in Dev**

- **Variant reporting** using custom scripts

**Note** Some of the later functions i.e. variant annotation and qc reporting are still in dev.  

**We highly recommed read trimming prior to alignment.**
We have noticed considerable speed-ups in alignmnet time and increased quality of SNP/INDEL calls using trimmed vs raw fastq.  

**Base quality score recalibration is also recommended.**  
As an alternative to GATK, we have added fucntionality for use of
**[BamUtil](https://github.com/statgen/bamUtil):[recab](http://genome.sph.umich.edu/wiki/BamUtil:_recab)**
for base quality score recalibration.

**Non-GATK users**
- are encouraged to use aligners such as **[stampy](http://www.well.ox.ac.uk/project-stampy)** and **[novoalign](http://www.novocraft.com)** that perform base quality score recal on the fly.  
- are encouraged to use variant callers that perform local re-aligmnet around candidate sites to mitigate the need for the indel realignment stages.  
    - **[freebayes](https://github.com/ekg/freebayes)**
    - **[platypus](http://www.well.ox.ac.uk/platypus)**

**************

## Dockerised NGS Tools

All NGSeasy Docker images can be pulled down from **[compbio Docker Hub](https://hub.docker.com/u/compbio/)** or using the Makefile.  
We provide an Amazon EBS data volume with indexed genomes: XXXXXX  


## Table 1. NGSeasy Tools

| Docker Image | Version | NGS Tool (version) | Short Description | URL |
|--------------|---------|--------------------|-------------------|-----|
| compbio/ngseasy-base | 1.0-r001 | VCFtools (v0.1.12b) | manipulate vcf | [link](http://vcftools.sourceforge.net/) |
| - | - | vt (latest) |  manipulate vcf | [link](http://genome.sph.umich.edu/wiki/Vt) |
| - | - | bcftools (1.2-5-g7fa0d25) |  manipulate vcf | [link](https://github.com/samtools/bcftools) |
| - | - | vcflib (v1.0.0) |  manipulate vcf | [link](https://github.com/ekg/vcflib) |
| - | - | samtools (1.2-17-ge91985a) |  manipulate sam/bam | [link](https://github.com/samtools/) |
| - | - | samblaster (0.1.21) |  manipulate sam/bam | [link](https://github.com/GregoryFaust/samblaster) |
| - | - | sambamba (v0.5.1)|  manipulate sam/bam | [link](https://github.com/lomereiter/sambamba) |
| - | - | bamUtil (1.0.13)|  manipulate sam/bam | [link](https://github.com/statgen/bamUtil) |
| - | - | bedtools (v2.23.0-10-g447cb97)|  manipulate bed files| [link](https://github.com/arq5x/bedtools2) |
| - | - | seqtk (1.0-r77-dirty)|  manipulate fastq| [link](https://github.com/lh3/seqtk) |
| - | - | vawk (0.0.2)|  manipulate vcf| [link](https://github.com/cc2qe/vawk) |
| - | - | bioawk (latest)|  manipulate sam/bam/vcf | [link](https://github.com/lh3/bioawk) |
| compbio/ngseasy-fastqc | 1.0-r001 | fastqc (v0.11.2) | FASTQ Quality Control Plots| [link](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) |
| compbio/ngseasy-trimmomatic | 1.0-r001 | trimmomatic (0.32) | FASTQ Quality Trimming | [link](http://www.usadellab.org/cms/index.php?page=trimmomatic) |
| compbio/ngseasy-bwa | 1.0-r001 | bwa ( 0.7.12-r1039) | Aligner | [link](http://bio-bwa.sourceforge.net/) |
| compbio/ngseasy-stampy | 1.0-r001 | stampy (stampy-1.0.27) | Aligner | [link](http://www.well.ox.ac.uk/project-stampy) |
| compbio/ngseasy-snap | 1.0-r001 | snap-aligner (1.0beta.18) | Aligner | [link](http://snap.cs.berkeley.edu/) |
| compbio/ngseasy-bowtie2 | 1.0-r001 | bowtie2 (2.2.4) | Aligner | [link](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) |
| compbio/ngseasy-novoalign | 1.0-r001 | novoalign (3.02.13) | Aligner | [link](http://www.novocraft.com/products/novoalign/) |
| compbio/ngseasy-gatk | 1.0-r001 | gatk (3.4-0) | NGS PowerTools | [link](https://www.broadinstitute.org/gatk/) |
| compbio/ngseasy-picardtools | 1.0-r001 | picardtools (1.128) | NGS PowerTools | [link](http://broadinstitute.github.io/picard/) |
| compbio/ngseasy-glia | 1.0-r001| glia (latest) | NGS local realignment | [link](https://github.com/ekg/glia) |
| compbio/ngseasy-platypus | 1.0-r001 | platypus (0.8.1) | Variant Caller | [link](https://github.com/andyrimmer/Platypus) |
| compbio/ngseasy-freebayes | 1.0-r001 | freebayes (v0.9.21-19-gc003c1e) | Variant Caller | [link](https://github.com/ekg/freebayes) |

## to add
- ABRA

## Running an NGSeasy Tool Interactively

Run as non-root user `pipeman`.

`-v /media/Data:/home/pipeman` : Mounts local directory `/media/Data` to container directory `/home/pipeman`

```bash
TOOL="bwa"

docker run \
-P \
-w /home/pipeman \
-e HOME=/home/pipeman \
-e USER=pipeman \
--user pipeman \
-v /media/Data:/home/pipeman \
-it compbio/ngseasy-${TOOL}:1.0 /bin/bash
```

************

Dockerised NGSeasy
==========================

![docker](https://github.com/KHP-Informatics/ngs/blob/master/figs/Docker_container_engine_logo.png)  

The following section describes getting the Dockerised NGSeasy Pipeline(s) and Resources, project set up and running NGSeasy.

Getting all resources and building required tools will take a few hours depending on network connections and any random "ghosts in the machine" - half a day in reality.
But once you're set up, thats it - you are good to go.

## System Requirements

See Table [**System Requirements**]() for our recommended system requirements.NGSeasy will run on any modern computer/workstation or cloud infrastructure.  The Hard Disk requirements are based on our experience and result from the fact that the pipeline/tools produce a range of intermediary and temporary files for each sample.

The full NGSeasy install includes indexed genomes for hg19 and b37 for all aligners, annotation files from GATK resource, and all of the NGSeasy docker images. Additional disk space is needed if the user wishes to install the databases associated with the variant annotators, Annovar, VEP and snpEff.

Based on our experience, a functional basic NGS compute system for a small lab, would consist of at least 4TB disk space, 60GB RAM and at least 32 CPU cores. Internet speed and network connectivity are a major bottle neck when dealing with NGS sized data, and groups are encouraged to think about these issues before embarking on multi sample or population level studies - where compute requirements can very quickly escalate.  

**System Requirements**

Component | Minimum | Recommended
|----|----|----|
RAM | 16GB | 48-60GB
CPU |  8 cores | 16-36 cores
Hard Disk (per sample) | 50-100GB | 200-500GB
NGSeasy Install | 200GB | 500GB
Annotation Databases | 500GB | >1TB


## Installing Docker

Follow the simple instructions in the links provided below  

- [Mac](https://docs.docker.com/installation/mac/)  
- [Windows](https://docs.docker.com/installation/windows/)
- [Ubuntu](https://docs.docker.com/installation/ubuntulinux/)

A full set of instructions for multiple operating systems are available on the [Docker website](https://docs.docker.com/installation/).

## Getting NGSeasy

We provide a simple **Makefile** to pull all of the public nsgeasy components, scripts and set up to correct project directory structre on your local machines.

Setting up the initial project can take up a day, depending on your local network connections and speeds.

The default install dir is the users **${HOME}** directory.
The  **Makefile** provides options to install to any user defined directory and select NGSeasy version. eg :-

```bash
## EG. Installing to /media/scratch
make INSTALLDIR="/media/scratch" VERSION="1.0" all
```

The  **Makefile** also allows installation of selected components (check out its insides!).


## Set up NGSeasy Project configuration file

Using Excel or something, make a **[config.file.tsv]**  file and save as [TAB] a Delimited file with ``.tsv`` extenstion.
This sets up Information related to: Project Name, Sample Name, Library Type, Pipeline to call, NCPU.

We provide a template that can be used with NGSeasy, see [ngseasy_test.config.tsv](https://docs.google.com/spreadsheets/d/13DosazuGdeAojZQ6YySM420p76RGleCYjR7a_MvcP2U/edit?usp=sharing).

The **[config.file.tsv]** should contain the following 23 columns for each sample to be run through a pipeline:-

|Variable|type|Description|Options(Examples)|
|--------|--------|--------|--------|
PROJECT_ID|STRING|Project ID| Cancer
SAMPLE_ID|STRING|Sample ID| SAMPLE_I
FASTQ1|STRING|Read 1 Fastq| foo_R1.fq.gz
FASTQ2|STRING|Read 2 Fastq| foo_R2.fq.gz
PROJECT_DIR|STRING|ngseasy project dir| /media/scratch/ngs_projects
DNA_PREP_LIBRARY_ID|STRING|NGS Library|
NGS_PLATFORM|STRING|NGS Platform|ILLUMINA
NGS_TYPE|STRING|NGS Type|WEX (exome), WGS (genome), TGS (targeted)
BAIT|STRING|bait bed file| FOO.bed
CAPTURE|STRING|Capture bed file| BAR.bed
GENOMEBUILD|STRING|genome verison|hg19, b37 , b38 (coming soon)
FASTQC|STRING|Select fastqc|no-fastqc, qc-fastqc
TRIM|STRING|Select trimming|no-trimm, atrimm, btrimm
BSQR|STRING|Select BSQR| no-bsqr, bam-bsqr, gatk-bsqr
REALN|STRING|Select Realignment| no-realn, bam-realn, gatk-realn
ALIGNER|STRING|Select Aligner|no-aln, bwa, stampy, snap, novoalign, bowtie2
VARCALLER|STRING|Select Variant Caller|no-varcall, freebayes, platypus, UnifiedGenotyper, HaplotypeCaller, ensemble
CNV|STRING|Select CNV caller|no-sv,all-sv,lumpy,delly,slope,exomedepth,mhmm,cnvnator
ANNOTATOR|STRING|Select variant annotator|no-anno,snpeff,annovar,vep
CLEANUP|STRING|clean up temp files|TRUE, FALSE
NCPU|NUMBER|number of cores|1 .. N
VERSION|NUMBER|NGSeasy version|1.0
NGSUSER|STRING|user email|stephen.j.newhouse@gmail.com

## The NGSeasy project directory
The user needs to make the relevent directory structures on their local machine before starting an NGS run.

On our sysetm we typically set up a top-level driectory called **ngs_projects** within which we store output from all our individual NGS projects.

Here we are working from local top level directory called **media/**, but this can really be any folder on your local system ie your home directory **~/${USER}**.  

Within this directory **media** we make the following folders: -

```bash
ngs_projects  
|  
|__raw_fastq  
|__config_files  
|__ngseasy_resources  
   |  
   |__reference_genomes_b37  
   |__reference_genomes_hg19
```

Running the script `make XXXX` ensures that all relevant directories are set up, and also enforces a clean structure to the NGS project.  

Within this we make a `raw_fastq` folder, where we temporarily store all the raw fastq files for each project. This folder acts as an initial stagging area for the raw fastq files. During the project set up, we copy/move project/sample related fastq files to their own specific directories.
Fastq files must have suffix and be gzipped: **\_1.fq.gz** or **\_2.fq.gz**  
furture version will allow any format  

Running `ngseasy` with the relevent configuration file, will set up the following directory structure for every project and sample within a project:-  

```bash
.
ngs_projects  
|  
|__raw_fastq  
|__config_files  
|__run_logs
|__ngseasy_resources
|
|__ project_id  
	|  
	|__run_logs  
	|__config_files  
	|
	|__sample_id_1  
	|	|  
	|	|__fastq  
	|	|__tmp  
	|	|__alignments  
	|	|__vcf  
	|	|__reports  
	|	|__config_files  
	|
	|
	|__sample_id_n  
		|  
		|__fastq  
		|__tmp  
		|__alignments  
		|__vcf  
		|__reports  
		|__config_files  
```

## The raw_fastq Directory

The ``raw_fastq`` Directory is a very special directory indeed.
This is where the user should copy and or move **ALL NEW RAW FASTQ** Files to.
This is to be used as an intial staging area for all fastq files.
NGSeasy expects all raw fastq data to be placed here for all new samples or runs.
NGSeasy inspects this folder and looks for the fastq file names specified in your confifg file.
If NGSeasy doen't find them, then it exits.
We do this to force the user to get organised.

****************

## Manually Build required NGSeasy Container Images

**Work In Progress...**

Currently we are not able to automatically build some of the tools in pre-built docker containers due to licensing restrictions.

Some of the software has restrictions on use particularly for commercial
purposes. Therefore if you wish to use this for commercial purposes, then you
leagally have to approach the owners of the various components yourself!  

**Software composing the pipeline requiring registration:-**  

   * novoalign http://www.novocraft.com/  
   * GATK https://www.broadinstitute.org/gatk/  
   * ANNOVAR http://www.openbioinformatics.org/annovar/  

**These tools require manual download and registration with the proivder. For non-academics/commercial groups, you will need to pay for some of these tools.**

Once you have paid/registered and downloaded the tool, we provide scripts and guidance for building these tools on your system.  

Its as easy as:-  

```bash
docker build -t compbio/ngseasy-${TOOL} .
```


## Building NOVOALIGN

**Download Novoalign from  http://www.novocraft.com/** into the local build directory **ngseasy/containerized/ngs_docker_debian/ngs_aligners/ngseasy_novoalign*.
Edit the [Dockerfile](https://github.com/KHP-Informatics/ngseasy/blob/master/containerized/ngs_docker_debian/ngs_aligners/ngseasy_novoalign/Dockerfile)
to relfect the correct version of novoalign.  

To use all novoalign fucntionality, you will need to **pay for a license**.

Once you obtained your **novoalign.lic**, download this to the build directory **ngseasy/containerized/ngs_docker_debian/ngs_aligners/ngseasy_novoalign*, which now should contain your updated [Dockerfile](https://github.com/KHP-Informatics/ngseasy/blob/master/containerized/ngs_docker_debian/ngs_aligners/ngseasy_novoalign/Dockerfile).

```bash
# move to ngseasy_stampy folder
cd ngseasy/containerized/ngs_docker_debian/ngs_aligners/ngseasy_novoalign
ls
```

**the directory should contain the following:-**

```
Dockerfile
novoalign.lic
README.md
novosortV1.03.01.Linux3.0.tar.gz
novocraftV3.02.08.Linux3.0.tar.gz
```

**build novoalign**

```bash
# build
docker build -t compbio/ngseasy-novoalign:v1.0 .
```

******

## Building GATK

You need to register and accept the GATK license agreement at https://www.broadinstitute.org/gatk/.  

Once done, download GATK and place in the GTAK build directory **ngseasy/containerized/ngs_docker_debian/ngs_utils/ngseasy_gatk**.  

Edit the [Dockerfile](https://github.com/KHP-Informatics/ngseasy/blob/master/containerized/ngs_docker_debian/ngs_utils/ngseasy_gatk/Dockerfile) to relfect
the correct version of GATK.  

```bash
# move to ngseasy_gatk folder
cd ngseasy/containerized/ngs_docker_debian/ngs_utils/ngseasy_gatk
ls
```

**the directory should contain the following:-**

```
Dockerfile
README.md
GenomeAnalysisTK-3.3-0.tar.bz2
```

**build gatk**

```bash
# build
docker build -t compbio/ngseasy-gatk:v1.0 .
```

********

## Manually Build NGSeasy Variant Annotaion Container Images

The tools used for variant annotation use large databases and the docker images exceed 10GB. Therefore, the user should manually build these container images prior to running the NGS pipelines.
Docker build files ([Dockerfile](https://docs.docker.com/jsearch/?q=Dockerfile)) are available for
- [Annovar](https://github.com/KHP-Informatics/ngs/tree/master/containerized/ngs_docker_debian/ngseasy_annovar/Dockerfile)  
- [VEP](https://github.com/KHP-Informatics/ngs/tree/master/containerized/ngs_docker_debian/ngseasy_vep/Dockerfile)
- [snpEff](https://github.com/KHP-Informatics/ngs/tree/master/containerized/ngs_docker_debian/ngseasy_snpeff/Dockerfile)  

**Note** Annovar requires user registration.  

**Once built on the user system, these container images can persist for as long as the user wants.**  

**Large Variant Annotation Container Images**

Its as easy as:-  

```bash
docker build -t compbio/ngseasy-${TOOL} .
```
********
## Build VEP
```bash

cd /media/ngs_projects/nsgeasy/ngs/containerized/ngs_docker_debian/ngseasy_vep

sudo docker build -t compbio/ngseasy-vep:${VERSION} .
```
********
## Build Annovar

```bash
cd /media/ngs_projects/nsgeasy/ngs/containerized/ngs_docker_debian/ngseasy_annovar

sudo docker build -t compbio/ngseasy-annovar:${VERSION} .
```
********
## Build snpEff
```bash
cd /media/ngs_projects/nsgeasy/ngs/containerized/ngs_docker_debian/ngseasy_snpeff

sudo docker build -t compbio/ngseasy-snpeff:${VERSION} .
```

********
## Coming Soon
- New Aligners:- [SNAP](http://snap.cs.berkeley.edu/), GSNAP, mr- and mrs-Fast,gem
- https://github.com/amplab/snap
- [SLOPE (CNV fo targetted NSG)] ((http://www.biomedcentral.com/1471-2164/12/184))
- Cancer Pipelines
- Annotation Pipelines and Databases
- Visualisation Pipelines
- Var Callers:- VarScan2
- SGE scripts and basic BASH scrips for running outside of Docker
- biobambam https://github.com/gt1/biobambam  
- bamaddrg https://github.com/ekg/bamaddrg  
- bamtools https://github.com/ekg/bamtools  

## Useful Links

- https://bcbio.wordpress.com/  
- https://basecallbio.wordpress.com/2013/04/23/base-quality-score-rebinning/  
- https://github.com/statgen/bamUtil  
- http://genome.sph.umich.edu/wiki/BamUtil:_recab  
- https://github.com/chapmanb/bcbio.variation  
- http://plagnol-lab.blogspot.co.uk/2013/11/faq-and-clarifications-for-exomedepth-r.html

*************************
(C) 2015 Kings College London

******

<p align="center">
  Development funded as part of:</br>
  <a href="http://www.maudsleybrc.nihr.ac.uk/">NIHR Maudsley Biomedical Research Centre (BRC), King's College London</a> and the </br>  
    <a href="http://www.ucl.ac.uk/health-informatics/">Farr Institute of Health Informatics Research, UCL Institute of Health Informatics, University College London</a>
  </br>
  <a href="http://phidatalab.org/">PHI Data Lab</a>:</br>
  <a href="https://www.kcl.ac.uk/ioppn/index.aspx">Institute of Psychiatry, Psychology & Neuroscience,King's College London</a>.
</p>
