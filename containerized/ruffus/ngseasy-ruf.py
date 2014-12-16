from ruffus import *

#########################################################################
# -- Author: Amos Folarin                                               #
# -- Author: Stephen J Newhouse                                         #
# -- Author: Kerz Maximilian                                            #
# -- Organisation: KCL/SLaM                                             #
# -- Email: amosfolarin@gmail.com                                       #
# -- Email: stephen.j.newhouse@gmail.com                                #
# -- Email: kerzmaximilian@gmail.com                                    #
#########################################################################


#------------------------------------------------------------------------
# A ruffus pipeline to manage the NGSeasy containerized pipeline
# TODO:
# 1) test general docker calling from ruffus
# 2) create a lib for ruffus pipeline function(s) which should map to one
#    container --> this will constitute a module
# 3) Build Web front end, drag & drop interface for construction of 
#    pipelines from these containerized modules
#------------------------------------------------------------------------



#------------------------------------------------------------------------
# pre-process/QC => {[fastqc] + [fastqc_d ruffus function]}
#------------------------------------------------------------------------




#------------------------------------------------------------------------
# bowtie module => {[bowtie container] + [bowtie_d ruffus function]}
#------------------------------------------------------------------------

@transform(starting_files               # Input = starting files
                    suffix(".fastq"),   # input file suffix = .fastq
                    ".bam")             # output file suffix = .bam

def bowtie_d(input_file,                 # 1st parameter is Input
                    output_file,        # 2nd parameter is Output
                    extra_params):      # user supplied args
    """
       Call bowtie docker container:  
          docker run
                    -v path/to/data/:/ngs/path/to/data
                    -v path/to/output:/ngs/path/to/output
                    ngseasy-bowtie bowtie [options] [extra_params]
                    input_file, output_file

        * although in reality I will probably wrap in a sepate func.
    """
    ii = open(input_file)
    oo = open(output_file, "w")
    oo.write(ii.read())



#------------------------------------------------------------------------
#
#------------------------------------------------------------------------

#------------------------------------------------------------------------
#
#------------------------------------------------------------------------

#------------------------------------------------------------------------
#
#------------------------------------------------------------------------

#------------------------------------------------------------------------
#
#------------------------------------------------------------------------

#------------------------------------------------------------------------
#
#------------------------------------------------------------------------

#------------------------------------------------------------------------
#
#------------------------------------------------------------------------












