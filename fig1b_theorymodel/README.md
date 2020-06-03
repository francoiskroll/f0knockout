# test
Data and scripts of preprint Kroll et al., 2020 – F0 knockout method

# about Illumina MiSeq analysis

AmpliCan requires fastq files as input, but I found that filtering my reads prior so I only use high-quality reads really improved the results.
As it is easier to filter reads on the .bam aligment file (eg. based on

## 1. Alignment + Filtering + Conversion

Reads were received as one forward fastq file & one reverse fastq file for each sample.
Script `alignFilterBackSingle.command` handles alignment/filtering/conversion back to fastq files, for one sample at a time (i.e. one pair of input fastq files).

(I use MacOS Terminal.)

So you can run the Shell (.command) script:

    chmod u+x ~/.../alignFilterBackSingle.command

Then, for example for sample

    alignFilterBack
