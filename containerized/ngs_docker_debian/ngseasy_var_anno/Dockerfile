# Base image
FROM compbioi/ngseasy-base:wheezy

# Maintainer 
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# Set correct environment variables.
ENV HOME /root

# Update
RUN apt-get update -y


#---------------------------------------------annotation-------------------------------------------------
# + SnpEff
RUN wget -O /tmp/snpEff_latest_core.zip http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip \
    && unzip /tmp/snpEff_latest_core.zip -d /usr/local/pipeline/ \
    && chown -R pipeman:ngsgroup /usr/local/pipeline/snpEff \
    && sed -i '$aCLASSPATH=.:${CLASSPATH}:/usr/local/pipeline/snpEff/snpEff.jar' /home/pipeman/.bashrc \
    && sed -i '$aPATH=${PATH}:/usr/local/pipeline/snpEff' /home/pipeman/.bashrc \
    && echo "alias ngsSNPeff='/usr/local/pipeline/snpEff'" >> /home/pipeman/.bashrc

# + ANNOVAR (see licence, and registration)
# Available on reg:  http://www.openbioinformatics.org/annovar/download/Ht8qRwQSTi/annovar.latest.tar.gz
RUN wget -O /tmp/annovar.latest.tar.gz http://www.openbioinformatics.org/annovar/download/mP628pfL21/annovar.latest.tar.gz \
  && tar xzvf /tmp/annovar.latest.tar.gz -C /usr/local/pipeline/ \
  && cd /usr/local/pipeline/annovar \
  && cp -v /usr/local/pipeline/annovar/table_annovar.pl table_annovar.pl.original \
  && wget http://www.openbioinformatics.org/annovar/download/table_annovar.pl \
  && chown -R pipeman:ngsgroup /usr/local/pipeline/annovar \
  && sed -i '$aPATH=${PATH}:/usr/local/pipeline/annovar' /home/pipeman/.bashrc \
  && echo "alias ngsSNPeff='/usr/local/pipeline/annovar'" >> /home/pipeman/.bashrc
  
#-------------------------------PERMISSIONS--------------------------
RUN chmod -R 766 /usr/local/pipeline/***
RUN chown -R pipeman:ngsgroup /usr/local/pipeline

# Cleanup the temp dir
RUN rm -rf /tmp/*

# open ports private only
EXPOSE 80

# Clean up APT when done.
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

