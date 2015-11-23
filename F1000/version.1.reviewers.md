# Brad Chapman
Bioinformatics Core, Department of Biostatistics, Harvard T.H. Chan School of Public Health, MA, USA

- Approved with Reservations

The authors describe NGSeasy, a containerized set of tools for variant calling on high throughput sequencing data. The implementation is open-source, maintained and well documented with instructions for getting started. The paper reports a number of useful practical results, including reproducibility through containerization, validations and timing of analyses. Below are my suggestions for helping to improve the paper and questions about the implementation details.

## Paper content
In the Introduction, you should note that similar Docker-based approaches exist for making Galaxy installation easier: https://github.com/bgruening/docker-galaxy-stable

You mention running at non-root in the Docker container but should also discuss that users need root privileges to install and run Docker. This currently limits usability of Docker on shared computing environments since it requires giving NGSeasy root-equivalent permissions. User namespace support is in progress and will help, but is not yet in any released versions of Docker. I agree with your security points once you have the NGSeasy approach setup, but getting there can be a challenge. You do mention this later in "NGSeasy future developments" and that would fit better in the initial installation section.

For reporting of download/install times, please also list install times from more standard connection speeds. A majority of users will not have 500Mb/s or better download. Is it possible to download subsets of the data? It looks like it currently grabs both hg19, b37 and hs37d5, tripling the download times and space required. Digging into the code it wasn't clear how to get other mentioned genomes like hs38DH.
Validation

For the GCAT/Genome in a Bottle validations, I'd suggest reporting precision instead of specificity. Specificity is not especially useful for calling since it's dominated by true negatives. For example, the precision rates show clear differences between FreeBayes and Platypus, and also differences between novoalign and the other aligners. The specificity numbers do not reveal these.

It's hard to judge the results of your validation without comparing to another best-practice pipeline like bwa + GATK HaplotypeCaller. Having these as a baseline next to your comparisons would strengthen the argument that the current implementation does a comparable job to expected best practice.

It would be useful to have bwa-mem alignment results also listed in the GCAT validations. bwa-mem is a widely used aligner, separate from stampy.

Do you have validation of using non-GATK tools (recab and glia) versus GATK tools in terms of the output quality? This would be useful to report. I've had good output success avoiding these step entirely but would like to see differences between avoiding the steps and using freely available alternatives.

## Timings
The timing information is really useful and a great addition to the paper. I'd suggest adding some caveats to the conclusion and tables to make it clearer about the inputs, since the numbers are exome with only 30x coverage. Most standard exomes would be higher coverage and WGS is becoming increasingly standard. Some of the statements like "alignment and variant calling are no longer a major bottle neck" seem overextended from timings on this smaller dataset. Scaling up is not linear and things get harder for WGS projects like 100k genomes project.

Has there been scaling work across non-single machine setups? Our experience is that shared network issues and managing Docker containers can dominate scaling. If the target is single multi-core machines it would be worth specifying this directly.
Docker
What is your experience with larger Docker containers and Docker Hub? Practically I've found a lot of timeout issues trying to download and manage larger images. Do you have workaround/experience with these issues?

Have you successfully run workflows on non-unix systems with Boot2Docker? You list these as workstations where NGSeasy should work. We've not had good success with mounting external data into Boot2Docker instances but would have interest in re-exploring if this changed: https://github.com/chapmanb/bcbio-nextgen-vm#mac-osx-docker-support

## Minor

Typo: sudo make INTSALLDIR="/media/scratch" all

I have read this submission. I believe that I have an appropriate level of expertise to confirm that it is of an acceptable scientific standard, however I have significant reservations, as outlined above.
Competing Interests: I work on the bcbio community project (https://github.com/chapmanb/bcbio-nextgen) which has overlapping aims to NGSeasy. I'm hopeful to collaborate on future NGSeasy development.

*************


# Fabien Campagne
Department of Physiology and Biophysics, Weill Cornell Medical College, NY, USA  

- Approved with Reservations

This manuscript describes a software tool called NGSEasy, which consists of a set of configured docker images for writing pipelines to call genomic variations. The manuscript addresses a need because the portability and reproducibility of bioinformatics pipelines, such as the one described in this manuscript, are very real problems faced by bioinformaticians and users of these methods. As such, the manuscript is an effort to describe an approach that could help with these challenges. Despite my interest in the manuscript, I have several reservations regarding the presentation of the results that I believe should be addressed before the manuscript is accepted for publication in a peer-reviewed journal. These reservations are:

The abstract indicates that the manuscript “demonstrates” best practices for packaging and executing multi-component pipelines for NGS. I am very uneasy about the term “best practices” used in a scientific context. I believe that it is very difficult to define objectively what constitutes “best practices”. Asking different experts may yield different answers, and one could argue that “best practices” are a synonym for asking an expert about his or her opinion. I believe that the manuscript would be strengthened if it described the authors’ recommendations and supported each recommendation with a clear and detailed rationale (perhaps outlining alternatives that the authors have tried and eventually rejected while developing NGSEasy, and explaining the reasons to do so).

While the word “demonstrates” in the abstract suggested to me that the manuscript would describe a set of practices recommended by the authors for packaging docker containers, it appears from the Result section that no such practices are explicitly described. It is therefore possible that the authors meant to refer to a set of published practices for analysis of high-throughput sequence data. If it is the case, as suggested by the sentence at the top right of page 4, I am not sure as to what is claimed in the abstract: demonstration of published practices, or recommended practices for packaging NGS code in a docker image?

The abstract and introduction define the domain of application of NGSEasy as “next-generation sequencing (NGS)”. However, the manuscript is about methods for variant calling, which is an important, but smaller scope that the full NGS data analysis. For instance, NGSEasy does not include tools for analysis of RNA-Seq. I recommend to revise the abstract and introduction to clearly indicate the scope of the software tool.

A strength of the manuscript is to use the GCAT server to evaluate the pipeline, but the results are not presented in the context of the performance of other pipelines, so the readers have no easy way of knowing if the sensitivity and specificity measures presented on page 6 are competitive. For instance, on Page 6, the manuscript claims: “we have successfully Dockerized a full NGS pipleline that is capable of producing meaningful result, that are comparable with public and “best practice” workflows”. However, there is no reference for the workflows the work is compared to and no simple way to establish if the results are comparable, let alone competitive. I strongly recommend to include a comparison directly in the manuscript to help the readers objectively assess performance.

The manuscript would be strengthened by providing a discussion of the limitations of the work. For instance, it is unclear what support is provided for parallelization across nodes, rather than SMP parallelization. (Multi-node parallelization is important when more than one or two samples need to be analyzed.)

I am unable to locate Reference 37 using the citation information:  “37. Matzke M, Jurkschat K, Backhaus T, et al.: PrePrints PrePrints. 2014; (1): 1–34. ”. This reference is used when discussing performance of docker containers and I am unable to determine if this is appropriate. A valid reference for this point is https://peerj.com/articles/1273/.

The reference provided for Nextflow is wrong. The tool should be cited using the web site (http://nextflow.io) or FigShare poster, and the correct authors given credit.  

## Minor comments:
Page 3, “NGSEasy contains all the basic tools needed for manipulation and quality control..” should be toned down. Using all in a manuscript is inviting contradiction. For instance, I could point out that the NGSEasy do not contain SpeedSeq, a recently published set of tools that considerably accelerates variation calling in HTS data. Therefore, I would argue that NGSEasy does not contain all the basic tools that I would like to use. Consider revising as “NGSEasy contains a set of tools sufficient for manipulation and quality control..”

Page 5. The word “all” is used again (left column, 6th paragraph). I doubt that the practice, as described, eliminates all potential issues with typo, since end-users will be writing scripts using tools in the image, and I am not sure how consistent naming conventions can fully eliminate typos when writing scripts.

Page 6. last paragraph, last sentence, check the grammar (missing a “y”?).

I have read this submission. I believe that I have an appropriate level of expertise to confirm that it is of an acceptable scientific standard, however I have significant reservations, as outlined above.
Competing Interests: No competing interests were disclosed

*************

# Michael Barton
 Joint Genome Institute, CA, USA

- Approved with Reservations

My understanding of this article is that NGSeasy pipeline aims to simplify the distribution of common tools used in sequencing analysis. A still significant problem in bioinformatics is getting the third-party tools installed and working, by using Docker containers as described in this article, the authors will make this process easier. The code is available on github as described and they provide extensive documentation.

## Major
One concern is the install instructions in the article. Specifically:

sudo make INTSALLDIR="/media/scratch" all
sudo make intsall

I am wary of using 'sudo' to install. I know that using tools like 'apt-get' require 'sudo' however for most bioinformatics software I prefer to install in my user directory simply to avoid any possible security problems. I took a look at the Makefile and I believe that sudo is not necessarily required to install, only the the INSTALLDIR and TARGET_BIN are owned by the user. Also there is typo here in 'intsall'

The project doesn't include any NGS tools related to assembly or transcriptomics. Though not stated specifically, the tools and data described here leads me to believe this project is focused around clinical applications and human genomics. If that is the case perhaps this should be clarified in the article and the title.

The Docker security issue at the end feels tagged-on. This, I think, is a pressing concern that prevents many people from using Docker on HPC machines as opposed to on-demand computing such as AWS. This is the case at the JGI where I currently work. I would suggest expanding on this point a little more to describe why it is an issue.

I think a short paragraph would be useful to end the article with. This would summarise the points described above and potential impact of the work.

I think Figure 1 could be expanded upon. It currently assumes a familiarity with Docker that could make it difficult to interpret without a good understanding of containers and volumes.

## Minor
The scripts installed from https://github.com/KHP-Informatics/ngseasy/tree/master/bin

The scripts here are very large (400-800 line) bash scripts. My concern is these scripts may be difficult to maintain or debug. Writing from my own experience, when projects get larger it is worth considering if code can be refactored or made more modular to make it more maintainable.

Inconsistent use of Dockerized / Dockerised between the article and the documentation

I believe Docker have deprecated boot2docker and now recommend docker-toolkit. This is their recommended way of installing Docker.

All data is available from the AWS EU region. This may take much longer to download outside of this region. I'm not sure if the authors can do anything about this however as using a CDN may be prohibitive.

The URL https://registry.hub.docker.com/repos/compbio requires the creation of a DockerHub account to view. This is not the developers fault, however having to register is likely to result in some users not following up.

I have read this submission. I believe that I have an appropriate level of expertise to confirm that it is of an acceptable scientific standard, however I have significant reservations, as outlined above.

Competing Interests: The authors invited me to present at their Bio in Docker conference in London this yea
