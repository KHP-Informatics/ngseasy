
```bash
## get ngs_resource mapping
myresources=`dirname ${PROJECT_DIR}`
NGSResources="${myresources}/ngs_resources"

${DOCKER_RUN} \
-v ${PROJECT_DIR}:/home/pipeman/ngs_projects \
-v ${NGSResources}:/home/pipeman/ngs_resources \
```

```bash

make INSTALLDIR="/media/Data" ngsprojectdir
make INSTALLDIR="/media/Data" b37
make INSTALLDIR="/media/Data" hg19
make INSTALLDIR="/media/Data" dockerimages
sudo make INSTALLDIR="/media/Data" install
```