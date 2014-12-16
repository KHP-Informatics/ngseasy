#!/bin/bash

#ANNOVAR DATA BASES
/usr/local/pipeline/annovar/annotate_variation.pl --buildver hg19 --downdb seq /home/pipeman/ngs_projects/humandb/hg19_seq

/usr/local/pipeline/annovar/annotate_variation.pl --buildver hg19 --downdb --webfrom annovar refGene  /home/pipeman/ngs_projects/humandb/

/usr/local/pipeline/annovar/annotate_variation.pl --buildver hg19 --downdb --webfrom annovar knownGene  /home/pipeman/ngs_projects/humandb/

/usr/local/pipeline/annovar/annotate_variation.pl --buildver hg19 --downdb --webfrom annovar ensGene  /home/pipeman/ngs_projects/humandb/

/usr/local/pipeline/annovar/retrieve_seq_from_fasta.pl /home/pipeman/ngs_projects/humandb/hg19_refGene.txt -seqdir /home/pipeman/ngs_projects/humandb/hg19_seq -format refGene -outfile /home/pipeman/ngs_projects/humandb/hg19_refGeneMrna.fa

/usr/local/pipeline/annovar/retrieve_seq_from_fasta.pl /home/pipeman/ngs_projects/humandb/hg19_knownGene.txt -seqdir /home/pipeman/ngs_projects/humandb/hg19_seq -format knownGene -outfile /home/pipeman/ngs_projects/humandb/hg19_knownGeneMrna.fa

/usr/local/pipeline/annovar/retrieve_seq_from_fasta.pl /home/pipeman/ngs_projects/humandb/hg19_ensGene.txt -seqdir /home/pipeman/ngs_projects/humandb/hg19_seq -format ensGene -outfile /home/pipeman/ngs_projects/humandb/hg19_ensGeneMrna.fa
