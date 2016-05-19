FROM snewhouse/ngseasybase:aplha-0.0.2

MAINTAINER Stephen Newhouse <stephen.j.newhouse@gmail.com>

# update and get some things
RUN apt-get update --fix-missing && \
    apt-get upgrade -y && \
    apt-get autoremove -y && \
    apt-get autoclean && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

## install ngs tools
ADD ngs_conda_tool_list.txt ngs_conda_tool_list.txt

RUN conda install -y --file ngs_conda_tool_list.txt

RUN nextflow self-update

CMD [ "/bin/bash" ]
