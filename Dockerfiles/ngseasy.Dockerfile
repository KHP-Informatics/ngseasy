FROM snewhouse/ngseasybase:aplha-0.0.5

MAINTAINER Stephen Newhouse <stephen.j.newhouse@gmail.com>

ADD ngseasy_conda_install.sh ngseasy_conda_install.sh

RUN /bin/bash ngseasy_conda_install.sh && \
rm ngseasy_conda_install.sh

VOLUME /home/ngseasy/anaconda2/bin

CMD [ "/bin/bash" ]
