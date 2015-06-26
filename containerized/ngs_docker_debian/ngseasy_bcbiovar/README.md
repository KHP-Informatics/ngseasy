# bcbio.variation

A toolkit to analyze genome variation data, built on top of the
[Genome Analysis Toolkit (GATK)][1] with Clojure. It supports scoring for the
[Archon Genomics X PRIZE competition][5] and is also a general framework for
variant file comparison. It enables validation of variants and exploration of
algorithm differences between calling methods by automating the process involved
with comparing two sets of variants. For users, this integrates with the
[bcbio-nextgen][8] framework to automate variant calling and validation. For
developers, bcbio.variation provides command line tools and an API to clean and
normalize variant data in [VCF format][2] avoiding comparison artifacts
associated with different variant representations.

* [Description of the comparison framework and variant calling algorithm comparisons][7]
* [Presentation from Bioinformatics Open Source Conference 2012][p1]
* [Presentation overview of the project][4]
* [Howto description of interfacing with GATK][6]
* [Code documentation][3]

[![Build Status](https://secure.travis-ci.org/chapmanb/bcbio.variation.png)](http://travis-ci.org/chapmanb/bcbio.variation)

[1]: http://www.broadinstitute.org/gsa/wiki/index.php/The_Genome_Analysis_Toolkit
[2]: http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40
[3]: http://chapmanb.github.com/bcbio.variation
[4]: http://chapmanb.github.com/bcbio.variation/presentations/gatk_clojure.pdf
[5]: http://genomics.xprize.org/
[6]: http://bcbio.wordpress.com/2012/03/04/extending-the-gatk-for-custom-variant-comparisons-using-clojure/
[p1]: http://chapmanb.github.com/bcbio.variation/presentations/variation_bosc_2012/variation_chapman.pdf
[7]: http://bcbio.wordpress.com/2013/05/06/framework-for-evaluating-variant-detection-methods-comparison-of-aligners-and-callers/
[8]: https://github.com/chapmanb/bcbio-nextgen

## Obtaining

### Download

The latest release is 0.1.9 (6 November 2014): [bcbio.variation-0.1.9-standalone.jar][dl].
Run from the command line:

    $ java -jar bcbio.variation-VERSION-standalone.jar [arguments]

The jar contains a full GATK commandline with additional walkers, as well as
custom command line programs. The usage section below contains examples of using
the library for variant comparison, normalization and ensemble calling.  Note
that bcbio.variation requires Java 1.7 since the underlying GATK libraries are
not compatible with earlier versions.

[dl]: https://github.com/chapmanb/bcbio.variation/releases/download/v0.1.9/bcbio.variation-0.1.9-standalone.jar
 
### As a library

To use as a library from Leiningen or Maven, follow the instructions on the
[clojars page][clojars].
 
[clojars]: https://clojars.org/bcbio.variation

### Latest source

The latest version is available directly from GitHub and requires Java 1.6 or
better and [Leiningen][u1] (version 2). Lein will automatically pull in all
required dependencies.

## Usage

### Generate summary of concordance between variant calls

To compare two GATK compatible VCF files in a specific region:

    $ java bcbio.variation.jar variant-utils comparetwo eval.vcf refcall.vcf ref.fa regions.bed

You can also tune many parameters with a [YAML][u2] configuration file specifies
the variant files for comparison. The project contains example configuration and
associated variant files that demonstrate the features of the library and the
configuration below has a description of available options.

An example of [scoring a phased diploid genome against a haploid reference genome][grade-ref]:

    $ java bcbio.variation.jar variant-compare config/reference-grading.yaml

An example of [assessing variant calls produced by different calling algorithms][grade-alg]:

    $ java bcbio.variation.jar variant-compare config/method-comparison.yaml

[grade-ref]: https://github.com/chapmanb/bcbio.variation/blob/master/config/reference-grading.yaml
[grade-alg]: https://github.com/chapmanb/bcbio.variation/blob/master/config/method-comparison.yaml

### Normalize a variant file

A tricky part of variant comparisons is that VCF format is flexible enough to
allow multiple representations. As a result two files may contain the same
variants, but one might have it present in a multi-nucleotide polymorphism while
another represents it as an individual variant.

To produce a stable, decomposed variant file for comparison run:

    $ java bcbio.variation.jar variant-prep your_variants.vcf your_reference.fasta

This will also handle re-ordering variants to match the reference file ordering,
essential for feeding into tools like GATK, and remapping hg19 to GRCh37
chromosome names. To additionally filter outputs by indel size, pass an argument
specifying the maximum indel size to include: `--max-indel 30`. To retain
reference (`0/0`) and no calls in the prepped file, use `--keep-ref`.

### Ensemble variant calls

bcbio.variation contains approaches to merge variant calls from multiple
approaches. It normalizes all input VCFs, prepares a combined variant file with
variants from all approaches, and then filter the combined file to produce a
final set of calls.

To combine multiple variant calls into a single combined ensemble callset:

    $ java bcbio.variation.jar variant-ensemble params.yaml reference.fa
                               out.vcf in1.vcf in2.vcf in3.vcf

- `params.yaml` -- parameters to us in preparing the ensemble set. See
  [the example config][e1] for options.
- `reference.fa` -- The reference FASTA file.
- `out.vcf` -- Name of the final combined dataset.
- `in1.vcf`... -- The input variant callsets to combine.

[e1]: https://github.com/chapmanb/bcbio.variation/blob/master/config/ensemble/combine-callers.yaml

### Web interface

The [o8 visualization framework][w1] provides a web interface to this toolkit,
providing interaction with Galaxy and GenomeSpace, visualization of biological
metrics associated with variants, reactive filtering and automated scoring.

[w1]: https://github.com/chapmanb/o8

### Run GATK walker for variant statistics

    $ lein uberjar
    $ java -jar target/bcbio.variation-0.0.1-SNAPSHOT-standalone.jar -T VcfSimpleStatsWalker
      -R test/data/GRCh37.fa --variant test/data/gatk-calls.vcf --out test.png

### Run custom GATK annotator

    $ lein uberjar
    $ java -jar target/bcbio.variation-0.0.1-SNAPSHOT-standalone.jar -T VariantAnnotator
       -A MeanNeighboringBaseQuality -R test/data/GRCh37.fa -I test/data/aligned-reads.bam
       --variant test/data/gatk-calls.vcf -o annotated-file.vcf

## Configuration file

A YAML configuration file defines targets for comparison processing. Two example
files for [reference grading][u4] and [comparison of calling methods][u3]
provide example starting points and details on available options are below:

    dir:
      base: Base directory to allow use of relative paths (optional).
      out: Working directory to write output.
      prep: Prep directory where files will be pre-processed.
    experiments: # one or more experiments
     - sample: Name of current sample.
       ref: Reference genome in FASTA format.
       intervals: Intervals to process in BED format (optional).
       align: Alignments for all calls in BAM format (optional).
       summary-level: Amount of summary information to provide,
                      [full,quick] (default:full)
       params: # Processing parameters associated with this experiment
         max-indel: Maximum indel size to include in non-structural variant
                    comparisons (default: 30)
         multiple-thresh: Threshold for percentage of comparisons to match
                          to consider two multiple sample variants the same.
                          (default: 1.0)
         compare-approach: Provide alternative approaches to compare variants:
                           approximate -- allow flexible matching of het/hom variants
                                          and indels
         recall-approach: Method to use for recalling variants:
                          consensus -- Most common variant from multiple inputs
       approach: Type of comparison to do [compare,grade]. Default compare.
       calls: # two or more calls to compare
         - name: Name of call type
           file: One or more input files in VCF format
           align: Alignment for specific call in BAM format (optional).
           ref: Reference genome if different than experiment ref (optional)
           intervals: Genome intervals to process in BED format (optional).
           metadata: Dictionary of annotations associated with the call set.
                     Finalizers use these to provide annotation specific
                     filtering of calls.
           filters: Provide hard filtering of variants prior to comparison with 
                    specified JEXL GATK expressions.
           format-filters: Provide hard filtering of variants based on
                           attributes in the genotype FORMAT field.
           recall: Recall variant positions after merging multiple input calls. 
                   Can tune with `recall-approach` in `params` (boolean; default false)
           annotate: Annotate calls with GATK annotations (boolean; default false).
           normalize: Normalize MNPs and indels (boolean: default true).
           prep: Prep with in-order chromosomes and sample names (boolean; default false).
           prep-sort-pos: Sort by position during prep. Required if variants are
                          not coordinate sorted within chromosomes. (boolean; default false).
           fix-sample-header: Adjust VCF sample header names to match sample
                              specified in `sample` (boolean; default false)
           prep-sv-genotype: Normalize structural variant genotypes to a single
                             ref call (boolean; default false).
           prep-allele-count: Number of alleles to convert calls to during
                              prep work (default 2)
           preclean: Remove problematic characters from input VCFs
                     (boolean; default false). 
           remove-refcalls: Remove reference, non-variant calls.
                            (boolean; default false). 
           make-haploid: Convert a set of diploid calls to haploid variants
                        (boolean; default false)

## Finalizer configuration

In addition to the pairwise comparisons, the configuration allows specification
of additional filtration and all-by-all comparisons based on the pairwise
results. Like `calls`, specify these under an experiment with the `finalize`
tag. Available methods are:

* `multiple` which does a comparison of a target call method to all other
  calls. A comparison of GATK calls to all other methods looks like:

        finalize:
          - method: multiple
            target: gatk
            ignore: []

  and produces three output files:

  - true positives -- calls detected in all methods
  - false negatives -- calls not found in gatk, but detected in all other methods
  - false positives -- calls found in gatk but callable and discordant
   in one of the other methods

  The `ignore` option provides a list of methods to ignore in the all-by-all
  overlap comparison.

* `recal-filter` to do post-comparisons filtering of calls. This can use either
  the results of a pairwise comparison or `multiple` comparison. An example
  demonstrating all of the filtering options re-filters a GATK versus FreeBayes
  comparison:

        finalize: 
          - method: recal-filter
            target: [gatk, freebayes]
            params:
              filters: [HRun > 5.0]
              classifiers: [AD, DP, QUAL]
              xspecific: true
              trusted:
                total: 0.75
                technology: 0.65
              untrusted:
                total: 0.25
    
The options for filtering are:

  - `filters` -- Perform hard filtering of the file with specified expressions.
  - `classifiers` -- Perform classification of true/false reads
   based on the provided attributes.
  - `trusted` -- Metadata annotation values that specify trusted variants not
   subjected to filtering. The example retains variants present in more than 75%
   of calls or 65% of different technologies.
  - `untrusted` -- Specify threshold for variants that should be automatically
   filtered. The example excludes variants with support from less than 25% of the
   callers.
  - `xspecific` -- Identify and filter calls specific to a technology or calling
    method, when combining multiple callers and technologies.

  You can specify the background to use for training with `support`. There are
  two options:

  - `support: gatk` -- Use an all-by-all comparison based on GATK to establish
   true and false positives.
  - `support: [gatk, freebayes]` -- Use the gatk/freebayes pairwise comparison
   for true and false positives.
           
[u1]: https://github.com/technomancy/leiningen#installation
[u2]: http://en.wikipedia.org/wiki/YAML
[u3]: https://github.com/chapmanb/bcbio.variation/blob/master/config/method-comparison.yaml
[u4]: https://github.com/chapmanb/bcbio.variation/blob/master/config/reference-grading.yaml
[u5]: http://www.broadinstitute.org/gsa/wiki/index.php/Variant_quality_score_recalibration

## Utilities

This library also contains useful command line utilities to help with variant
preparation and analysis:

- Create a BAM file compatible with GATK. This converts coordinates between hg19
  and GRCh37 for human samples if needed, reorders chromosomes relative to the
  input file and adds run group information with a defined sample name:

        $ lein variant-reorder your_file.bam /path/to/GRCh37.fa SampleName
    
- Provide a summary CSV file of call information for a VCF file, including
  mappings back to an original set of pairwise analyses:

        $ lein variant-utils callsummary variants.vcf original-combined-config.yaml

- Sort a variant VCF file to reference ordering defined in a FASTA file

        $ lein variant-utils sort-vcf variants.vcf reference.fa

- Convert an Illumina directory of variant calls into a single, cleaned VCF:

        $ lein variant-utils illumina /path/to/IlluminaDir sample-name GRCh37.fa hg19.fa

## Contributors

- Brad Chapman
- Chris Fields
- Kevin Lynagh
- Justin Zook

## License

The code is freely available under the [MIT license][l1].

[l1]: http://www.opensource.org/licenses/mit-license.html