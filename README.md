NGSeasy-v1.0 (beta/dev)
===================
**A [Dockerized](https://www.docker.com/) and [Virtulaized](https://www.virtualbox.org/) ngs pipeline and tool-box.** 

**With NGSeasy you can now have full suite of NGS tools up and running on any high end workstation in an afternoon**

**Note:** NGSeasy is under continous development and the dev version evolves quickly. NGSeasy-v1.0 Full Production release will be available Dec 2014


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
- Indexed Reference Genomes
- Cancer Pipelines
- Annotation Pipelines
- Visualisation Pipelines

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

![ngsEASY](figs/ngsEASY_atomic_pipeline_visualisation.png "Dockerized NGS Pipeline")

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


Running the Dockerised NGSeasy Pipeline
==========================================
**With NGSeasy and the simple steps outlined below, you can now have the full suite of NGS tools up and running on any high end workstation in an afternoon**

The commands below walk you through getting, setting up and running an NGSeasy pipeline.  

**Example Set up and  Running NGSeasy**

**1.** Make a directory on your local machine for storing data, NGSeasy scripts and run files generated from the pipeline.

```sh
# make ngseasy directory
mkdir /media/ngseasy
```

**2.** Get the latest NGSeasy pipeline scripts from GitHub

```sh
# get latest ngseasy_scripts code
cd /media/ngseasy
git clone https://github.com/KHP-Informatics/ngs.git
cd ngs
# move scripts to top level folder
mv -v ngseasy_scripts /media/ngseasy/
```

**3.** Get the latest NGSeasy Resources

```sh
cd  /media/ngseasy/

# get ngseasy_resources (on Dropbox for a limited time only) or sftp (as above)
wget --no-check-cert https://www.dropbox.com/s/9pw3ml75pdnufjl/ngseasy_resources.tar.gz

# Un pack resources 
tar -xrfv ngseasy_resources.tar.gz

# move [fastq_raw], [reference_genomes_b37] and [gatk_resources] to /media/ngseasy/
mv -v /media/ngseasy/ngseasy_resources/reference_genomes_b37 /media/ngseasy/
mv -v /media/ngseasy/ngseasy_resources/gatk_resources /media/ngseasy/
mv -v /media/ngseasy/ngseasy_resources/gatk_resources /media/ngseasy/
mv -v /media/ngseasy/ngseasy_resources/fastq_example /media/ngseasy/
```

**4.** Make ``ngs_projects`` NGS projects directory
```sh
# move up to ngseasy
cd /media/ngseasy

# make ngs_projects
mkdir /media/ngseasy/ngs_projects
```

**On your local machine, you should now have the following directories:-**

- ``fastq_raw`` [To hold all incoming raw fastq files]
- ``ngs_projects`` [Output directory for all NGS projects/runs]
- ``reference_genomes_b37`` [b37 Reference Genomes indexed for BWA/Bowtie2/Stampy/Novovalign]
- ``gatk_resources`` [vcf and annotation files used by GATK]
- ``ngseasy_scripts`` [copy of latest pipeline scripts]
- ``ngs`` [GitHub Clone]  

**NB:** These folders are required by pipeline as they are hardcoded in the NGSeasy scripts.

**5.** Copy ``run_ngseasy_dockers.sh`` to ``ngs_projects`

``run_ngseasy_dockers.sh`` calls and runs the NGSeasy pipelines

```sh
# copy run_ngseasy_dockers.sh to ngs_projects
cp -v /media/ngseasy/ngseasy_scripts/run_ngseasy_dockers.sh /media/ngseasy/ngs_projects
```

**6.** Set Up Config File and copy to ``ngs_projects``

```sh
cp -v ngseasy_scripts/ngs.config.file.tsv /media/ngseasy/ngs_projects
```

In Excel make config file and save as [TAB] Delimited file with ``.tsv`` extenstion.  
See Example provided and [GoogleDoc](https://docs.google.com/spreadsheets/d/1kp1Nyw0x3zXqO2Wm2Z25ErJ0Z-Uoab8tjRPq9h4sonk/edit?usp=sharing). Remove the header from this file before running the pipeline. This sets up Information related to: Project Name, Sample Name, Library Type, Pipeline to call, NCPU.

The [config.file] should contain the following 15 columns for each sample to be run through a pipeline:- 

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


**7.** Get Main NGSeasy Docker Image (if you haven't already done this!)

```sh
sudo docker pull compbio/ngseasy-alignment-public:v1.2
```

**8.** Run an NGSeasy Pipeline

```sh
cd /media/ngseasy

bash run_ngseasy_docker.sh -c ngs.config.file.tsv -d /media/ngseasy

```


[Back to The Begining](https://github.com/KHP-Informatics/ngs/blob/master/containerized/README.md#ngs-easy-v10)

Putting it all together
----------------------------

```bash
# Make NGSeasy Directory
$ mkdir /media/ngseasy

# Get NGSeasy scripts
$ cd /media/ngseasy
$ git clone https://github.com/KHP-Informatics/ngs.git
$ cp -v ./ngs/ngseasy_scripts /media/ngseasy

# get ngseasy_resources (on Dropbox for a limited time only) or sftp (as above)
$ wget --no-check-cert https://www.dropbox.com/s/9pw3ml75pdnufjl/ngseasy_resources.tar.gz?dl=0;
$ tar -xrvf ngseasy_resources.tar.gz
$ mv -v ./ngseasy_resources/** /media/ngseasy
# copies to [fastq_raw], [reference_genomes_b37] and [gatk_resources] to ngseasy/
# Make local directory for NGS Projects
$ mkdir /media/ngseasy/ngs_projects

# In Excel, setup config file and save as tab delimited text file *.tsv. See Example provided. This sets up Information related to: Project Name, Sample Name, Library Type, Pipeline to call, NCPU etc etc

# Get Main NGSeasy Docker Image
$ sudo docker pull compbio/ngseasy-alignment-public:v1.2

# Run NGSeasy
$ cd ngs_projects
$ bash run_ngseasy_docker.sh -c ngs.config.file.tsv -d /media/ngseasy

```

**And its as simple as that!**

[Back to The Begining](https://github.com/KHP-Informatics/ngs/blob/master/containerized/README.md#ngs-easy-v10)

**************************

NGSeasy Pipeline Scripts
=========================

**Main Pipeline Calling**
- run_ngseasy_dockers.sh
- run_ea-ngs.sh 

**Full Pipeline from Fastq to VCF**
- full_gatk.sh
- full_no_gatk.sh

**NGSeasy Components**
- fastqc_and_trimm.sh
- alignment.sh
- post_alignment_coverage.sh
- post_alignment_qc.sh
- variant_calling_multiple_samples.sh
- variant_calling_single_sample.sh
- sv-cnmops.sz
- sv-delly.sh
- sv-ensemble.sh
- sv-exomedepth.sh
- sv-lumpy.sh
- sv-mhmm.sh

[Back to The Begining](https://github.com/KHP-Informatics/ngs/blob/master/containerized/README.md#ngs-easy-v10)

***************************



************************

### Getting the GATK Resources Bundle for yourself

[What's in the resource bundle and how can I get it?](https://www.broadinstitute.org/gatk/guide/article.php?id=1213)

```bash
## FTP Login Details
location: ftp.broadinstitute.org
username: gsapubftp-anonymous
password: <blank>
```
I would recommend using a separate program like [FileZilla](https://filezilla-project.org/), which will make it much easier for you to set up and manage your file transfers

**Path to gatk_resources.tgz**
```
ftp://ftp.broadinstitute.org/distribution/gsa/gatk_resources.tgz
```

[Back to The Begining](https://github.com/KHP-Informatics/ngs/blob/master/containerized/README.md#ngs-easy-v10)

**********


NGSEASY-VM : An NGS Tool Box
================================
![VirtualBox](figs/add-virtualbox.png "VirtualBox")  

A virtual machine based on [Ubuntu 14.04 LTS](http://www.ubuntu.com/desktop), containing all these [NGS Tools](https://github.com/KHP-Informatics/ngs/blob/master/containerized/README.md#the-tools-included-are-as-follows--) (and a few extras) is available. The virtual machine contains indexed reference genomes (b37) for all the installed aligners and the [GATK Resources Budle](https://www.broadinstitute.org/gatk/download), and includes the latest version of Docker.

Using [VirtualBox](https://www.virtualbox.org/) and our NGSEASY-VM [ngseasy-vm-v1.0.vdi]() you can have a full suite of NGS tools up and running on any high end workstation in less than an hour.

## 1. Installing Oracle VirtualBox

The latest Oracle VM VirtualBox is available [here](http://www.oracle.com/technetwork/server-storage/virtualbox/downloads/index.html#vbox)

Remember to install the latest [Oracle VM VirtualBox Extension Pack](http://www.oracle.com/technetwork/server-storage/virtualbox/downloads/index.html#extpack)

A full set of instructions are available here:-  

- [Oracle VirtualBox](http://www.oracle.com/technetwork/server-storage/virtualbox/downloads/index.html)  
- [VirtualBox User Manual](http://download.virtualbox.org/virtualbox/UserManual.pdf)  
- [VirtualBox Installation Details](https://www.virtualbox.org/manual/ch02.html)  


## 2. Getting the ngseasy-vm

Download the ngseasy-vm [ngseasy-vm-v1.0.vdi]() and save this anywhere on your local machine.  

**NOTE:** This file is XXXGB in size.  

## 3. Installing the ngseasy-vm

Create a new virtual machine from the downloaded [ngseasy-vm-v1.0.vdi]().  

**Step 1.** Fire up VirtualBox and click **New**  
**Step 2.** The **Create Virtual Machine** window will appear  
**Step 3.** You will presented with 3 options that need to be set as  

> Name: **ngseasy**  
> Type: **Linux**  
> Version: **Ubuntu (64 bit)** 
 
**Step 4.** Click **Next** . This will take you to the **Memory Size** window  
**Step 5.** Select the amount of memmory (RAM) in megabytes to be allocated to the virtual machine (recommended 8GB). The size you 
allocate will depend on the available RAM on your machine. Adjust this vaule as needed. Click **Next**  
**Step 6.** This will take you to the **Hard drive** window. This is the important step. You will be presented with 3 options:-  

> Do not add a virtual hard drive  
> Create a virtual hard drive  
> Use an existing virtual hard drive file  

**Step 7.** Select **Use an existing virtual hard drive file**.   
**Step 8.** Click on the folder icon and point to the **ngseasy-vm.vdi** file.  
**Step 9.** Step  Click **Create**  

Once the new virtual machine is created, that virtual machine should fire up as expected 
and regardless of platform (I have tested this going from Linux to Linux, Linux to Windows, and Windows to Linux hosts).  

Before starting the virtual machine you can go to the **Settings** tab and alter the amount or RAM, CPUs assigned to the machine.
There is also the functionality to mount shared drives between the host OS and the virtual machine. In the **General** > **Advanced** tab,
ensure thet the **Shared Clipboard** and **Drang n Drop** are both set to **Bidirectional**.

[Back to The Begining](https://github.com/KHP-Informatics/ngs/blob/master/containerized/README.md#ngs-easy-v10)

******

Minimum Hardware Requirements
===============================
- 64 bit computer and operating system is recommended for more than 2 GB RAM
- 8GB RAM
- 8 Cores
- Minimum 100GB storage

[Back to The Begining](https://github.com/KHP-Informatics/ngs/blob/master/containerized/README.md#ngs-easy-v10)

******

```
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
```
[Back to The Begining](https://github.com/KHP-Informatics/ngs/blob/master/containerized/README.md#ngs-easy-v10)

******

