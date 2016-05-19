FROM snewhouse/ngseasybase:aplha-0.0.1

MAINTAINER Stephen Newhouse <stephen.j.newhouse@gmail.com>

ENV LANG en_US.UTF-8
ENV LC_ALL en_US.UTF-8

ADD ngs_conda_tool_list.txt ngs_conda_tool_list.txt
ADD requirements.txt requirements.txt

## update conda
RUN conda update conda
RUN conda install -y --file requirements.txt
RUN conda update conda-build 

## add channels
RUN conda config --add channels bioconda
RUN conda config --add channels r
RUN conda config --add channels sjnewhouse

## install ngs tools
RUN conda install -y --file ngs_conda_tool_list.txt

RUN nextflow self-update

CMD [ "/bin/bash" ]
