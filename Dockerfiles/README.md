## Dockerfiles for ngseasy images


### basimage
sets up ngseasy user and some of the volumes  
installs some basic tools to build stuff  

### ngseasy
- Uses [`conda`](http://conda.pydata.org/docs/) to add ngs tools (see `ngs_conda_tool_list.txt`)  
- Anaconda (https://www.continuum.io/downloads)  
- See install script for details. see `ngseasy_conda_install.sh`  
- The single "APP" here is an NGS tool box for alignment and variant calling  
- See `ngseasy_conda.Dockerfile`  
- Full Anaconda install
- miniconda maybe : after we work out all the required dependencies
- i have a feeling it will still be quite large image

**ngseasy_conda.Dockerfile**  

- test build  

```bash
docker build --force-rm --rm=true --file=ngseasy_conda.Dockerfile -t snewhouse/ngseasy:aplha-0.0.3 .
```

- v0.0.3 : after fixing mess with volumes being ROOT and adding anaconda2/bin to PATH
- **4.934 GB** Image

```
REPOSITORY          TAG                 IMAGE ID            CREATED             SIZE
snewhouse/ngseasy   aplha-0.0.3         47bde4210d45        5 seconds ago       4.934 GB
```

**docker-squash**

- https://github.com/jwilder/docker-squash

```bash
docker save 47bde4210d45 | sudo docker-squash -t snewhouse/ngseasy:squash-a-0.0.3 | docker load
```

## Create NGSeasy env?
This file `ngseasy-spec-file-22-05-16.txt` may be used to create an environment using:

```bash
conda create --name ngseasy --file ngseasy-spec-file-22-05-16.txt
```

## A thought

```bash
docker run \
-v ${DATA}:/home/ngseasy/ngs_projects \
-v ${RAWFASTQ}:/home/ngseasy/ngs_projects/fastq \
-v ${REFGENOMES}:/home/ngseasy/reference_genomes \
-v ${NOVOALIGN_LIC}:/home/ngseasy/anaconda2/bin \
-v ${GATK}:/home/ngseasy/anaconda2/bin \
--rm=true \
-it snewhouse/ngseasy:aplha-0.0.2 bash ngseasy -c runconfig.tsv
```

## Warning
conda builds for some important tools are missing for osx_64 Mac and I guess Windows. 

**Use Linux Folks** and run NGS on HPC or the cloud or Workstation with >16CPUs, > 30GB RAM and LOOOOTS of storage!

## misc

- dev on rosalind openstack takes toooo long.
- network issues mean that pulling images takes ages
- wget install file for anaconda > 24hours
- moved dev to AWS; using my personal account (t2.small , 50GB Storage)
- all built in less than an hour!!!!!!

### smaller images
- see mulled project: https://github.com/mulled/mulled
- future work
- need missing packages up on conda or bioconda first
- time!
