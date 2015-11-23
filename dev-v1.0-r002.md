## NGSeasy Re-Dev (F1000 revisions)

|Details|
|:----|
|Name: Dr Stephen Newhouse|
|Email: stephen.j.newhouse@gmail.com|
|Date:30.10.2015|
|Git: https://github.com/KHP-Informatics/ngseasy|

## Check List

- [x] Build base debain image: fat dev box image
- [x] Build Big Kahuna : All the tools in one place
- [ ] Build individual tools from `<tool>.Dockerfile`
  - [ ] base image small
  - [ ] fastqc
  - [ ] trimmomatic
  - [ ] picardtools
  - [ ] bwa
  - [ ] freebayes
  - [ ] gatk
  - [ ] snap
  - [ ] novoalign
  - [ ] bowtie2
  - [ ] stampy
  - [ ] mrsfast
  - [ ] platypus
  - [ ] samtools
  - [ ] VarScan
  - [ ] VarDict
- [ ] Build tool sets: Base, Aligners, Callers...
- [ ] Build speedseq
- [ ] Build bcbio-nextgen
- [ ] Get all genomes and decoys
 - [ ] hg19
 - [ ] b37
 - [ ] b38
 - [ ] Add decoy genomes to all  
 - [ ] get resources and lift-over all
- [ ] Index all genomes for all aligners
 - [ ] bwa
 - [ ] novo
 - [ ] snap
 - [ ] bowtie2
 - [ ] mrsfast  
- [ ] Test Suite for GCAT and Paper
- [ ] Add HLA Calling
- [ ] Add Gemini
- [ ] Add Nextflow Versions
- [ ] Add Native non-Docker version


## Dev Box:

```
Linux sjn-devbox-ubuntu 3.19.0-31-generic #36-Ubuntu SMP Wed Oct 7 15:04:02 UTC 2015 x86_64 x86_64 x86_64 GNU/Linux
```

### lscpu

```
ubuntu@sjn-devbox-ubuntu:~/sjn_bin$ lscpu
Architecture:          x86_64
CPU op-mode(s):        32-bit, 64-bit
Byte Order:            Little Endian
CPU(s):                32
On-line CPU(s) list:   0-31
Thread(s) per core:    1
Core(s) per socket:    1
Socket(s):             32
NUMA node(s):          1
Vendor ID:             GenuineIntel
CPU family:            6
Model:                 42
Model name:            Intel Xeon E312xx (Sandy Bridge)
Stepping:              1
CPU MHz:               2699.998
BogoMIPS:              5399.99
Hypervisor vendor:     KVM
Virtualization type:   full
L1d cache:             32K
L1i cache:             32K
L2 cache:              4096K
NUMA node0 CPU(s):     0-31
```

### Hardware

```
ubuntu@sjn-devbox-ubuntu:~/sjn_bin$ sudo lshw -short
H/W path        Device  Class      Description
==============================================
                        system     Computer
/0                      bus        Motherboard
/0/0                    memory     230GiB System memory
/0/1                    processor  Intel Xeon E312xx (Sandy Bridge)
/0/2                    processor  Intel Xeon E312xx (Sandy Bridge)
/0/3                    processor  Intel Xeon E312xx (Sandy Bridge)
/0/4                    processor  Intel Xeon E312xx (Sandy Bridge)
/0/5                    processor  Intel Xeon E312xx (Sandy Bridge)
/0/6                    processor  Intel Xeon E312xx (Sandy Bridge)
/0/7                    processor  Intel Xeon E312xx (Sandy Bridge)
/0/8                    processor  Intel Xeon E312xx (Sandy Bridge)
/0/9                    processor  Intel Xeon E312xx (Sandy Bridge)
/0/a                    processor  Intel Xeon E312xx (Sandy Bridge)
/0/b                    processor  Intel Xeon E312xx (Sandy Bridge)
/0/c                    processor  Intel Xeon E312xx (Sandy Bridge)
/0/d                    processor  Intel Xeon E312xx (Sandy Bridge)
/0/e                    processor  Intel Xeon E312xx (Sandy Bridge)
/0/f                    processor  Intel Xeon E312xx (Sandy Bridge)
/0/10                   processor  Intel Xeon E312xx (Sandy Bridge)
/0/11                   processor  Intel Xeon E312xx (Sandy Bridge)
/0/12                   processor  Intel Xeon E312xx (Sandy Bridge)
/0/13                   processor  Intel Xeon E312xx (Sandy Bridge)
/0/14                   processor  Intel Xeon E312xx (Sandy Bridge)
/0/15                   processor  Intel Xeon E312xx (Sandy Bridge)
/0/16                   processor  Intel Xeon E312xx (Sandy Bridge)
/0/17                   processor  Intel Xeon E312xx (Sandy Bridge)
/0/18                   processor  Intel Xeon E312xx (Sandy Bridge)
/0/19                   processor  Intel Xeon E312xx (Sandy Bridge)
/0/1a                   processor  Intel Xeon E312xx (Sandy Bridge)
/0/1b                   processor  Intel Xeon E312xx (Sandy Bridge)
/0/1c                   processor  Intel Xeon E312xx (Sandy Bridge)
/0/1d                   processor  Intel Xeon E312xx (Sandy Bridge)
/0/1e                   processor  Intel Xeon E312xx (Sandy Bridge)
/0/1f                   processor  Intel Xeon E312xx (Sandy Bridge)
/0/20                   processor  Intel Xeon E312xx (Sandy Bridge)
/0/100                  bridge     440FX - 82441FX PMC [Natoma]
/0/100/1                bridge     82371SB PIIX3 ISA [Natoma/Triton II]
/0/100/1.1              storage    82371SB PIIX3 IDE [Natoma/Triton II]
/0/100/1.2              bus        82371SB PIIX3 USB [Natoma/Triton II]
/0/100/1.2/1    usb1    bus        UHCI Host Controller
/0/100/1.2/1/1          input      QEMU USB Tablet
/0/100/1.3              bridge     82371AB/EB/MB PIIX4 ACPI
/0/100/2                display    GD 5446
/0/100/3        eth0    network    Virtio network device
/0/100/4                storage    Virtio block device
/0/100/5                generic    Virtio memory balloon
/0/100/6                storage    Virtio block device
/0/100/7                storage    Virtio block device
```


### Ubuntu version

```
LSB Version:	core-2.0-amd64:core-2.0-noarch:core-3.0-amd64:core-3.0-noarch:core-3.1-amd64:core-3.1-noarch:core-3.2-amd64:core-3.2-noarch:core-4.0-amd64:core-4.0-noarch:core-4.1-amd64:core-4.1-noarch:security-4.0-amd64:security-4.0-noarch:security-4.1-amd64:security-4.1-noarch
Distributor ID:	Ubuntu
Description:	Ubuntu 15.04
Release:	15.04
Codename:	vivid
```

## Docker Version

```
Client:
 Version:      1.9.0
 API version:  1.21
 Go version:   go1.4.2
 Git commit:   76d6bc9
 Built:        Tue Nov  3 17:48:04 UTC 2015
 OS/Arch:      linux/amd64

Server:
 Version:      1.9.0
 API version:  1.21
 Go version:   go1.4.2
 Git commit:   76d6bc9
 Built:        Tue Nov  3 17:48:04 UTC 2015
 OS/Arch:      linux/amd64
```

## java version

```bash
java version "1.7.0_85"
OpenJDK Runtime Environment (IcedTea 2.6.1) (7u85-2.6.1-6+deb8u1)
OpenJDK 64-Bit Server VM (build 24.85-b03, mixed mode)
```

## Python Version

```bash
Python 2.7.9
```

******************

## Aligners

```
aligners: bwa, snap, mrsfast, bowtie2, novoalign, stampy
```

| Aligners |
|:------|
|bwa|
|snap|
|mrsfast|
|bowtie2|
|novoalign|
|stampy|


## Basic Variant Callers
```
varcallers: freebayes, platypus, VarDict, VarScan, samtools,gatk-hc
```

|Variant Callers |
|:------|
| freebayes |
| platypus |
| VarDict |
| VarScan |
| samtools |
| gatk-hc |


## Structural Variant Callers

```
svcallers: lumpy,lumpyexpress,svtyper, cnvkit
```

|Structural Variant Callers |
|:------|
|lumpy |
|lumpyexpress|
|svtyper|
|cnvkit|

## Re-Aligners

```
realn: abra.jar, glia, ogap
```

|Re-Aligners|
|:------|
|abra|
|glia|
|ogap|

************

## To Do
- cwl  
- nextflow
- https://www.dockstore.org  
- Fat, Thin and Tiny Containers    

## Random Lists...
Some tools to add?  

- [List of sequence alignment software ](https://en.wikipedia.org/wiki/List_of_sequence_alignment_software)  
- http://bioinform.github.io/metasv/
