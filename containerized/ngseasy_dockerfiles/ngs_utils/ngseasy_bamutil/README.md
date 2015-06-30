# BamUtil and libStatGen

bamUtil is a repository that contains several programs that perform operations on SAM/BAM files. All of these programs are built into a single executable, bam.  

libStatGen  

### Full documentation 

http://genome.sph.umich.edu/wiki/BamUtil  
http://genome.sph.umich.edu/wiki/C%2B%2B_Library:_libStatGen  
### Git

https://github.com/statgen/bamUtil  


## Build Docker image

```bash
sudo docker build -t compbio/ngseasy-bamtuil .
```

## Get Docker image

```bash
sudo docker pull compbio/ngseasy-bamtuil:latest
```

## Run Docker image
- http://genome.sph.umich.edu/wiki/BamUtil:_recab  

```bash
sudo docker run -i -t compbio/ngseasy-bamtuil /usr/local/pipeline/bamUtil/bin/bam recab --in ${INPUT}.bam --out ${OUTPUT}.bam --refFile ${REF} --dbsnp ${DBSNP} --storeQualTag OQ --maxBaseQual 40 
```

