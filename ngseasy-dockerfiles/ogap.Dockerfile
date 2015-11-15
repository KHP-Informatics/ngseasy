

RUN cd /usr/local/ngs/bin && \
git clone --recursive https://github.com/ekg/ogap.git && \
cd ogap && \
make all && \
cp -v ogap /usr/local/bin

# ogap  and bamleftalign
 cd /usr/local/ngs/bin/ && \
    git clone --recursive https://github.com/ekg/ogap.git && \
    cd ogap && \
    make all && \
    chmod -R 777 ./* && \
    cp -v ogap /usr/local/bin/ && \
    cd /usr/local/ngs/bin/ && \
    git clone --recursive git://github.com/ekg/freebayes.git && \
    cd freebayes && \
    make all && \
    chmod -R 777 ./* && \
    cp bin/bamleftalign /usr/local/bin/ && \
    rm -frv /usr/local/ngs/bin/freebayes && \

USER ngseasy
WORKDIR /home/ngseasy
