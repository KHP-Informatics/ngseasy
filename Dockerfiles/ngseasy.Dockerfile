FROM snewhouse/ngseasybase:aplha-0.0.3

MAINTAINER Stephen Newhouse <stephen.j.newhouse@gmail.com>

## install ngs tools
ADD ngseasy_conda_local_install.sh ngseasy_conda_local_install.sh

RUN bash ngseasy_conda_local_install.sh && rm ngseasy_conda_local_install.sh

ENV PATH /home/ngseasy/anaconda2/bin

CMD [ "/bin/bash" ]
