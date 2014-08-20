NGS EASY v1.0
===================

[Dockerized](https://www.docker.com/) and [Virtulaized](https://www.virtualbox.org/) ngs pipeline and tool-box.  

As a multi-component system, NGS pipeline setup is traditionally heavy on 
configuration. Our idea is to provide this in a simple encapsulated container. 
Users also typically wish to configure their own environments and run the 
pipeline on a wide range of hardware (workstations to clusters to cloud), being 
able to stand-up a pipeline with minimal fuss is made straightforward with this 
container.  

**With NGSEASY you can now have full suite of NGS tools up and running on any high end workstation in less than an hour**

<a href="https://twitter.com/share" class="twitter-share-button">Tweet</a>
<script>!function(d,s,id){var js,fjs=d.getElementsByTagName(s)[0],p=/^http:/.test(d.location)?'http':'https';if(!d.getElementById(id)){js=d.createElement(s);js.id=id;js.src=p+'://platform.twitter.com/widgets.js';fjs.parentNode.insertBefore(js,fjs);}}(document, 'script', 'twitter-wjs');</script>

### Authors
- Amos Folarin <amosfolarin@gmail.com>  
<a href="http://www.linkedin.com/pub/amos-folarin/34/b06/978">
<img src="http://www.linkedin.com/img/webpromo/btn_viewmy_160x33.png" width="160" height="33" alt="View Amos's profile on LinkedIn">
</a>

- Stephen J Newhouse <stephen.j.newhouse@gmail.com>  
<a href="http://uk.linkedin.com/pub/dr-stephen-newhouse/29/89a/11a">
<img src="http://www.linkedin.com/img/webpromo/btn_viewmy_160x33.png" width="160" height="33" alt="View Steve's profile on LinkedIn">
</a>

- David Brawand <dbrawand@nhs.net>  
- Aditi Gulati <aditigulati3@googlemail.com>  

Lets us know if you want other tools added to NGSEASY  

## Table of Contents
[NOTICE TO USERS](https://github.com/KHP-Informatics/ngs/blob/dev/containerized/README.md#notice-to-users-of-the-container-image-or-vm)  

[Software requiring registration](https://github.com/KHP-Informatics/ngs/blob/dev/containerized/README.md#software-composing-the-pipeline-requiring-registration)  

[Overview of Pipeline Components](https://github.com/KHP-Informatics/ngs/blob/dev/containerized/README.md#overview-of-pipeline-components)  

[NGS Tools](https://github.com/KHP-Informatics/ngs/blob/dev/containerized/README.md#the-tools-included-are-as-follows--)  

[Dockerised NGSEASY set up](https://github.com/KHP-Informatics/ngs/blob/dev/containerized/README.md#dockerised-ngseasy)  

[Installing Docker](https://github.com/KHP-Informatics/ngs/blob/dev/containerized/README.md#installing-docker)  

[Getting the Dockerised NGSEASY Pipeline](https://github.com/KHP-Informatics/ngs/blob/dev/containerized/README.md#getting-the-dockerised-ngseasy-pipeline)  

[Running the Dockerised NGSEASY Pipeline](https://github.com/KHP-Informatics/ngs/blob/dev/containerized/README.md#running-the-dockerised-ngseasy-pipeline)  

[Local Machine Set up](https://github.com/KHP-Informatics/ngs/blob/dev/containerized/README.md#local-machine-set-up)  

[NGSEASY-VM : An NGS Tool Box](https://github.com/KHP-Informatics/ngs/blob/dev/containerized/README.md#ngseasy-vm--an-ngs-tool-box)  

[Installing Oracle VirtualBox](https://github.com/KHP-Informatics/ngs/blob/dev/containerized/README.md#installing-oracle-virtualbox)  

[Getting the ngseasy-vm](https://github.com/KHP-Informatics/ngs/blob/dev/containerized/README.md#getting-the-ngseasy-vm)  

[Installing the ngseasy-vm](https://github.com/KHP-Informatics/ngs/blob/dev/containerized/README.md#installing-the-ngseasy-vm)  

******

# NOTICE TO USERS OF THE CONTAINER IMAGE OR VM

While the software used to build the image is composed of free software versions
some of the software has restrictions on use particularly for commercial 
purposes. Therefore if you wish to use this for commercial purposes, then you 
leagally have to approach the owners of the various components yourself!  

This pipeline uses a number of pieces of software which require registration. 
By using this you are agreeing to observe the Terms and Conditions of the 
relevant pieces of software that compose this pipeline.  

# Software composing the pipeline requiring registration

If you want to build the image from the Dockerfile then you need to get your 
own versions of (below) in the build directory:

   * novoalign http://www.novocraft.com/
   * Stampy http://www.well.ox.ac.uk/project-stampy
   * GATK https://www.broadinstitute.org/gatk/
   * ANNOVAR http://www.openbioinformatics.org/annovar/

[Back to The Begining](https://github.com/KHP-Informatics/ngs/blob/dev/containerized/README.md#ngs-easy-v10)

******

Overview of Pipeline Components
================================
The basic pipeline contains all the basic tools needed for manipulation and 
quality control of raw fastq files (ILLUMINA focused), SAM/BAM manipulation,
alignment, cleaning (based on GATK best practises [ADD LINK]) and first pass
variant discovery. Separate containers are provided for indepth variant annotation,
structural variant calling, basic reporting and visualisations.  

![ngsEASY](figs/ngsEASY_pipeline_visualisation.png "Dockerized NGS Pipeline")


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
- lumpy-sv
- delly
- m-HMM
- CNVnator
- ExomeDepth

[Back to The Begining](https://github.com/KHP-Informatics/ngs/blob/dev/containerized/README.md#ngs-easy-v10)

******

##### To Add 
- SegSeq
- MuTect
- MutSig

******

Dockerised NGSEASY
==========================
![docker](figs/Docker_container_engine_logo.png "Docker")  

## Installing Docker

Follow the simple instructions in the links provided below  

- [Mac](https://docs.docker.com/installation/mac/)  
- [Windows](https://docs.docker.com/installation/windows/)
- [Ubuntu](https://docs.docker.com/installation/ubuntulinux/)

A full set of instructions for multiple operating systems are available on the [Docker website](https://docs.docker.com/installation/).

## Getting the Dockerised NGSEASY Pipeline

Available Dockerised NGSEASY containers at [Docker Hub](https://hub.docker.com/u/afolarin/)  
 
- afolarin/seq-alignment
- afolarin/var-calling
- afolarin/sv-calling
- afolarin/var-anno
- afolarin/vis-reports
- delly

```bash
sudo docker pull afolarin/seq-alignment
```

## Running the Dockerised NGSEASY Pipeline

```bash
sudo docker run \
-v ~/FASTQ_STAGGING:~/FASTQ_STAGGING \
-v ~/reference_geneomes:~/reference_genomes \
-v ~/ngs_projects:~/ngs_projects \
-u pipeman -t [CONTAINER] ngs.config
```

## Local Machine Set up

1. On your local machine, make the following directories:-

- FASTQ_STAGGING [FASTQ_STAGGING is to hold all incoming raw fastq files]
- ngs_projects [out put directory for all ngs projects)
- reference_genomes [get from URL and unpack. NB: XXX GB!]
- gatk_resources
- humandb [annovar]
- vep [VEP data base]
- snpeff [snpeff database]

These folders are required by pipeline as they are hardcoded in the nsg scripts.

- Get reference genomes, gatk resouces, snp annotation databases from [ADD URL]
- Un compress and save on local machine

[Back to The Begining](https://github.com/KHP-Informatics/ngs/blob/dev/containerized/README.md#ngs-easy-v10)

******

NGSEASY-VM : An NGS Tool Box
================================
![VirtualBox](figs/add-virtualbox.png "VirtualBox")  

A virtual machine based on [Ubuntu 14.04 LTS](http://www.ubuntu.com/desktop), containing all these [NGS Tools](https://github.com/KHP-Informatics/ngs/blob/dev/containerized/README.md#the-tools-included-are-as-follows--) (and a few extras) is available. The virtual machine contains indexed reference genomes (b37) for all the installed aligners and the [GATK Resources Budle](https://www.broadinstitute.org/gatk/download).  

Using [VirtualBox](https://www.virtualbox.org/) and our NGSEASY-VM you can have a full suite of NGS tools up and running on any high end workstation in less than an hour.

## 1. Installing Oracle VirtualBox

The latest Oracle VM VirtualBox is available [here](http://www.oracle.com/technetwork/server-storage/virtualbox/downloads/index.html#vbox)

Remember to install the latest [Oracle VM VirtualBox Extension Pack](http://www.oracle.com/technetwork/server-storage/virtualbox/downloads/index.html#extpack)

A full set of instructions are available here:-  

- [Oracle VirtualBox](http://www.oracle.com/technetwork/server-storage/virtualbox/downloads/index.html)  
- [VirtualBox User Manual](http://download.virtualbox.org/virtualbox/UserManual.pdf)  
- [VirtualBox Installation Details](https://www.virtualbox.org/manual/ch02.html)  


## 2. Getting the ngseasy-vm

Download the ngseasy-vm [ngseasy-vm.vdi]() and save this anywhere on your local machine.  

**NOTE:** This file is XXXGB in size.  

## 3. Installing the ngseasy-vm

Create a new virtual machine from the downloaded [ngseasy-vm.vdi]().  

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

[Back to The Begining](https://github.com/KHP-Informatics/ngs/blob/dev/containerized/README.md#ngs-easy-v10)

******

Minimum Hardware Requirements
===============================
- 64 bit computer and operating system is recommended for more than 2 GB RAM
- 8GB RAM
- 8 Cores
- Minimum 100GB storage

[Back to The Begining](https://github.com/KHP-Informatics/ngs/blob/dev/containerized/README.md#ngs-easy-v10)

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
[Back to The Begining](https://github.com/KHP-Informatics/ngs/blob/dev/containerized/README.md#ngs-easy-v10)

******

