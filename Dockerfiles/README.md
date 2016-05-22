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

**ngseasy_conda.Dockerfile**  

- test build  

```bash
docker build --force-rm --rm=true --file=ngseasy_conda.Dockerfile -t snewhouse/ngseasy:aplha-0.0.1 .
```

## Create NGSeasy env?
This file `ngseasy-spec-file-22-05-16.txt` may be used to create an environment using:

```bash
conda create --name ngseasy --file ngseasy-spec-file-22-05-16.txt
```

## Warning
conda builds for some important tools are missing for osx_64 Mac and I guess Windows. 

**Use Linux Folks** and run NGS on HPC or the cloud or Workstation with >16CPUs, > 30GB RAM and LOOOOTS of storage!

## misc

- dev on rosalind openstack takes toooo long.
- network issues mean that pulling images takes ages
- wget install file for anaconda > 24hours
- moved dev to AWS; using my personal account (t2.small , 50GB Storage)
