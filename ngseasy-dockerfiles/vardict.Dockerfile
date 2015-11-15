# base image
FROM compbio/ngseasy-base:r1.0-002
# Maintainer
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com
LABEL Description="This is the varDict images" Version="r1.0-002"

# varDict
RUN cd /usr/local/ngs/bin/ && \
    git clone --recursive https://github.com/AstraZeneca-NGS/VarDictJava.git && \
    cd VarDictJava && \
    git checkout v1.4.0 && \
    git submodule update --recursive && \
    cd /usr/local/ngs/bin/VarDictJava && \
    ./gradlew clean installApp && \
    chmod -R 755 /usr/local/ngs/bin/VarDictJava && \
    chown -R ngseasy:ngseasy /usr/local/ngs/bin && \
    sed  -i '$aPATH=$PATH:/usr/local/ngs/bin//VarDictJava/VarDict' /home/ngseasy/.bashrc
#./gradlew clean javadoc && \

CMD ["/bin/bash"]
