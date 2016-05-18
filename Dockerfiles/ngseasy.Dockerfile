FROM snewhouse/ngseasybase:aplha-0.0.1

MAINTAINER Stephen Newhouse <stephen.j.newhouse@gmail.com>

RUN mkdir -p /home/ngseasy/conda/conda-bld/linux-64 /home/ngseasy/conda/conda-bld/osx-64

ADD ngs_conda_tool_list.txt ngs_conda_tool_list.txt

RUN conda update -y conda
RUN conda update -y conda-build
RUN conda config --add channels bioconda
RUN conda config --add channels r
RUN conda config --add channels sjnewhouse

RUN conda install -y --file ngs_conda_tool_list.txt
RUN nextflow self-update

CMD [ "/bin/bash" ]
