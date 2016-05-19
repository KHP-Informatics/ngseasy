FROM snewhouse/ngseasybase:aplha-0.0.2

MAINTAINER Stephen Newhouse <stephen.j.newhouse@gmail.com>


RUN conda update -y conda
RUN conda update -y --all
RUN conda install -y \
  python=2.7 \
  pysam \
  pyyaml \
  java-jdk

## install ngs tools
ADD ngs_conda_tool_list.txt ngs_conda_tool_list.txt

RUN conda install -y \
--update-dependencies \
--file ngs_conda_tool_list.txt

RUN nextflow self-update

CMD [ "/bin/bash" ]
