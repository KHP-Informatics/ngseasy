## samtools

The original samtools package has been split into three separate
but tightly coordinated projects:
- [htslib](https://github.com/samtools/htslib): C-library for handling high-throughput sequencing data
- samtools: mpileup and other tools for handling SAM, BAM, CRAM
- [bcftools](https://github.com/samtools/bcftools): calling and other tools for handling VCF, BCF

See also http://github.com/samtools/

## HTSlib
HTSlib is an implementation of a unified C library for accessing common file
formats, such as [SAM, CRAM and VCF][1], used for high-throughput sequencing
data, and is the core library used by [samtools][2] and [bcftools][3].
HTSlib only depends on [zlib][4].
It is known to be compatible with gcc, g++ and clang.

HTSlib implements a generalized BAM index, with file extension `.csi`
(coordinate-sorted index). The HTSlib file reader first looks for the new index
and then for the old if the new index is absent.

This project also includes the popular tabix indexer, which indexes both `.tbi`
and `.csi` formats, and the bgzip compression utility.

[1]: http://samtools.github.io/hts-specs/
[2]: http://github.com/samtools/samtools
[3]: http://samtools.github.io/bcftools/
[4]: http://zlib.net/

## BCFtools
This is the official development repository for BCFtools. It contains all the vcf* commands
which previously lived in the htslib repository (such as vcfcheck, vcfmerge, vcfisec, etc.)
and the samtools BCF calling from bcftools subdirectory of samtools. 

For a full documentation, see [bcftools GitHub page](http://samtools.github.io/bcftools/). 

Other useful links:
------------------

File format specifications live on [HTS-spec GitHub page](http://samtools.github.io/hts-specs/)
[htslib](https://github.com/samtools/htslib)
[samtools](https://github.com/samtools/samtools)
[tabix](https://github.com/samtools/tabix)

