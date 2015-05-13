
```bash
## get ngs_resource mapping
myresources=`dirname ${PROJECT_DIR}`
NGSResources="${myresources}/ngs_resources"

${DOCKER_RUN} \
-v ${PROJECT_DIR}:/home/pipeman/ngs_projects \
-v ${NGSResources}:/home/pipeman/ngs_resources \
```