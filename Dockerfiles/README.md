## Dockerfiles for ngseasy images


### basimage
sets up ngseasy user and some of the volumes  
installs some basic tools to build stuff  

### ngseasy
- uses `conda` to add ngs tools    
- see install script for details  
- the single "APP" here is an NGS tool box for alignment and variant calling 

## Warning
conda builds for some important tools are missing for osx_64 Mac and I guess Windows. 

**Use Linux Folks** and run NGS on HPC or the cloud or Workstation with >16CPUs, > 30GB RAM and LOOOOTS of storage!

## misc

- dev on rosalind openstack takes toooo long.
- network issues mean that pulling images takes ages
- wget install file for anaconda > 24hours
- moved dev to AWS; using my personal account (t2.small , 50GB Storage)
