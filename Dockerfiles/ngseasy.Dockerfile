FROM snewhouse/ngseasybase:aplha-0.0.1

MAINTAINER Stephen Newhouse <stephen.j.newhouse@gmail.com>

RUN conda update -y conda && \
 conda update -y conda-build && \
 conda config --add channels r && \
 conda config --add channels bioconda && \
 conda config --add channels sjnewhouse

ADD ngs_conda_tool_list.txt /home/ngseasy

RUN conda install -y --file /home/ngseasy/ngs_conda_tool_list.txt && \
 nextflow self-update

CMD [ "/bin/bash" ]
