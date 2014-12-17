# ANNOVAR #

### ANNOVAR: Functional annotation of genetic variants from high-throughput sequencing data ###

ANNOVAR is an efficient software tool to utilize update-to-date information to functionally annotate genetic variants detected from diverse genomes (including human genome hg18, hg19, as well as mouse, worm, fly, yeast and many others). Given a list of variants with chromosome, start position, end position, reference nucleotide and observed nucleotides, ANNOVAR can perform:

Gene-based annotation: identify whether SNPs or CNVs cause protein coding changes and the amino acids that are affected. Users can flexibly use RefSeq genes, UCSC genes, ENSEMBL genes, GENCODE genes, or many other gene definition systems.
Region-based annotations: identify variants in specific genomic regions, for example, conserved regions among 44 species, predicted transcription factor binding sites, segmental duplication regions, GWAS hits, database of genomic variants, DNAse I hypersensitivity sites, ENCODE H3K4Me1/H3K4Me3/H3K27Ac/CTCF sites, ChIP-Seq peaks, RNA-Seq peaks, or many other annotations on genomic intervals.
Filter-based annotation: identify variants that are reported in dbSNP, or identify the subset of common SNPs (MAF>1%) in the 1000 Genome Project, or identify subset of non-synonymous SNPs with SIFT score>0.05, or find intergenic variants with GERP++ score<2, or many other annotations on specific mutations.
Other functionalities: Retrieve the nucleotide sequence in any user-specific genomic positions in batch, identify a candidate gene list for Mendelian diseases from exome data, and other utilities.
TABLE_ANNOVAR is a script within the ANNOVAR package that is very popular among users. Given a list of variants from whole-exome or whole-genome sequencing, it will generate an Excel-compatible file with gene annotation, amino acid change annotation, SIFT scores, PolyPhen scores, LRT scores, MutationTaster scores, PhyloP conservation scores, GERP++ conservation scores, dbSNP identifiers, 1000 Genomes Project allele frequencies, NHLBI-ESP 6500 exome project allele frequencies and other information.

In a modern desktop computer (3GHz Intel Xeon CPU, 8Gb memory), for 4.7 million variants, ANNOVAR requires ~4 minutes to perform gene-based functional annotation, or ~15 minutes to perform stepwise "variants reduction" procedure, making it practical to handle hundreds of human genomes in a day.

****

Thank you for your interests in using ANNOVAR. An updated version (2014nov12) is
available.  

Here is the link to download ANNOVAR:  
http://www.openbioinformatics.org/annovar/download/g4EUwyphi9/annovar.latest.tar.gz  

Please feel free to contact me if you have any questions or comments, but do NOT
reply this email as it will be IGNORED! Use a different title in your email.  

Please cite ANNOVAR if you use it in your research (Wang K, Li M, Hakonarson H.
ANNOVAR: Functional annotation of genetic variants from next-generation
sequencing data, Nucleic Acids Research, 38:e164, 2010). I spent tremendous
amount of time and efforts to maintain this tool, and your citation really means
a lot to me.  

### ADDITIONAL NOTES TO ANNOVAR USERS: ###  

Several notable improvements in the new version of ANNOVAR include: (1)
significantly reduce memory usage for filter annotation (2) improve
compatibility for unconventional chromosome names for species such as tomato (3)
fixed a problem in exon numbering for splice variants in reverse strand (4)
other minor functional improvements and database updates.  

Several notable additions/updates to ANNOVAR databases include: (1) updated
refGene, knownGene and ensGene files for hg18/hg19/hg38 are available to
download (2) 1000 Genomes Project 2014Oct version is available for all
populations and five subpopulations (AFR,AMR,EAS,EUR,SAS) with chrX/Y markers.
(3) ExAC 65000 exomes v0.2 is available all all populations and seven
subpopulations (AFR,AMR,EAS,FIN,NFE,OTH,SAS). (4) Updated ljb26 databases from
dbNSFP indexed by ANNOVAR is available to download. (5) Updated Clinvar and
Cosmic are available to download.  

Point of interest to ANNOVAR users: (1) A recent Genome Biology paper
(http://genomebiology.com/2014/15/3/R53) demonstrated that ANNOVAR is the most
widely used tool (>50% share) for interpretation of clinical genome sequencing
data (2) A recent Genome Medicine paper
(http://genomemedicine.com/content/6/3/26) reported that the choice of
transcript set can have a large effect on the ultimate variant annotations from
ANNOVAR and VEP (ENSEMBL’s Variant Effect Predictor). (3) Multiple users asked
my comments to a blog (http://blog.goldenhelix.com/?p=2486). I recently posted a
reply at the end of the page, clarifying that ANNOVAR does generate correct and
biologically meaningful results in all examples shown in the blog. (4) Some
users may not be aware of the wANNOVAR server (http://wannovar.usc.edu) which
allows annotation of VCF files via a web interface. Compared to ANNOVAR, it has
limited functionality and the data sources may not be as updated; however, it is
more user friendly.  

Mailing list related issues: (1) I typically send out an update email
every 4 months. (2) if you wish to be removed from the list, please email me
directly with "unsubscribe ANNOVAR" in the title of the email. (3) Starting from
next update, I may send emails to user via professional mailing list, rather
than from my personal email account.  

Finally, thank you very much for your help and support to ANNOVAR!  