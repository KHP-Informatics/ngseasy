NGSeasy NOTES
===============
- Stephen Newhouse <stephen.j.newouse@gmail.com>

*******************

**Get the Docker Imahge**
*****************************

```sh
sudo docker pull snewhouse/ngseasy-alignment-public:v1.2
```

**Get the pipeline scripts**
****************************

Get ``ngseasy_scripts``

```sh
git clone NNNNNNNNNNNNNN
```

Copy ``run_ngseasy_dockers.sh`` to ``ngs_projects``  
Copy ``config`` file to ``ngs_projects``  

**Set Up Config File**
************************

In Excel make config file and save as [TAB] Delimited file with ``.tsv`` extenstion.  
See Example.

**Call the pipeline**
********************
```sh
sh run_ngseasy_dockers.sh ngs.config.file.tsv
```

**Inside run_ngseasy_dockers.sh**
********************

```bash
#!/bin/bash
config_tsv=${1} 
sudo docker run -P \
-v /media/D/docker_ngs/ngseasy/fastq_raw:/home/pipeman/fastq_raw \
-v /media/D/docker_ngs/ngseasy/reference_genomes_b37:/home/pipeman/reference_genomes_b37 \
-v /media/D/docker_ngs/ngseasy/gatk_resources:/home/pipeman/gatk_resources \
-v /media/D/docker_ngs/ngseasy/ngs_projects:/home/pipeman/ngs_projects \
-v /media/D/docker_ngs/ngseasy/ngseasy_scripts:/home/pipeman/ngseasy_scripts \
-i -t snewhouse/ngseasy-alignment-public:v1.2 /sbin/my_init -- /bin/bash  /home/pipeman/ngseasy_scripts/run_ea-ngs.sh /home/pipeman/ngs_projects/${config_tsv};
```
Need **FULL PATHS** so ``../DIR`` should be ``/HOME/DIR``

### Misc



**Software Versions**
************************
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


**Trimmomatic Log**
********************
Specifying a trimlog file creates a log of all read trimmings, indicating the following details:

- the read name
- the surviving sequence length
- the location of the first surviving base, aka. the amount trimmed from the start
- the location of the last surviving base in the original read
- the amount trimmed from the end


```
206B4ABXX100825:7:1:7466:34224/1 69 0 69 7
206B4ABXX100825:7:1:7466:34224/2 62 0 62 14
206B4ABXX100825:7:1:7466:35353/1 69 0 69 7
206B4ABXX100825:7:1:7466:35353/2 76 0 76 0
206B4ABXX100825:7:1:7466:35395/1 54 0 54 22
206B4ABXX100825:7:1:7466:35395/2 76 0 76 0
206B4ABXX100825:7:1:7466:38189/1 52 0 52 24
206B4ABXX100825:7:1:7466:38189/2 76 0 76 0
206B4ABXX100825:7:1:7466:39218/1 76 0 76 0

```