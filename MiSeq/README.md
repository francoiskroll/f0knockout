# MiSeq

___

Please cite the preprint if you use some of the data or code! <br />
https://www.biorxiv.org/content/10.1101/2020.06.04.133462v3

___

Folders <br />
* **amplican_runs**
* **fastarefs**

can be found at the Zenodo version of this repository (too heavy for GitHub):

https://doi.org/10.5281/zenodo.3900611

## Files included

In *amplican_runs*, there is one folder per ampliCan run: *crystal* (i.e. *slc45a2* + *mitfa* + *mpv17*);  *csnk1db*; *off* (for *slc24a5* off-targets); etc.

In each run folder:

* **_README.txt**: Terminal commands to filter the fastq files prior to ampliCan
* **_config.csv**: config file used by ampliCan (as in: https://rdrr.io/bioc/amplican/f/vignettes/amplicanOverview.Rmd)
* **amplican_.R**: R script to run ampliCan
* **fastq** folder:
  * raw fastq files (*...fastq.gz*). Each sample has two fastq files: *R1* for forward reads, *R2* for reverse reads.
  * bam alignment files (name is the well from the original 96-well plates)
  * bam.bai index file
  * filtered bam file (*..._filter.bam*): bam file after filtering (see below)
  * index file for the filtered bam file (*..._filter.bam.bai*)
  * In folder **filterfastq**: fastq files after filtering, converted back from the filtered bam file (see below)

## Analysis until ampliCan

(I used MacOS Terminal)

ampliCan requires fastq files as input, but I found that filtering my reads prior analysis to keep only high-quality reads really improved the detection of mutations.
I found it easier to filter reads on the bam alignment file, so the idea is:

* 1– prepare the fasta reference for alignment;
* 2– align the fastq reads to the reference amplicon;
* 3– filter the bam alignment file;
* 4– convert the bam file back to one forward + one reverse fastq files
* 5– run ampliCan on the filtered fastq files

To run the *.command* scripts, you probably need to do first:

    chmod u+x ~/.../script.command

### 1– Prepare fasta reference

Reference amplicon sequence (5'–3') starts with the forward primer and ends with the reverse-complement of the reverse primer. All in lowercase, except PAM sequence in uppercase.

Reference sequences are in folder **fastarefs** (each file *ref_gene_locus.fa*), the other files are generated when preparing the reference to be used by bwa mem (not exactly sure which of these files are needed or not).

Script `prepareFastaRef.command` prepares the reference from each fasta file.

For example:

    prepareFastaRef.command ref_slc24a5_aa.fa

(the commands are included in each *README.txt*)

### 2/3/4– Align reads to reference + Filter alignment file + Convert back to filtered fastq files

MiSeq reads were received as one forward fastq file (*R1*) & one reverse fastq file (*R2*) for each sample.

Script `alignFilterBackSingle.command` handles all of alignment/filtering/conversion back to fastq files, for one sample at a time (i.e. one pair of input fastq files).

Usage is:

    alignFilterBackSingle.command FORWARD_READS REFERENCE

For example for sample A2 (*slc24a5*, locus AA):

        alignFilterBackSingle.command 20190805-A2_S2_L001_R1_001.fastq.gz ref_slc24a5_aa.fa

(the commands are included in each *README.txt*)

The script looks for the reverse reads by looking for a file of the same name but *R2* instead of *R1*.

Alignment algorithm is `bwa mem`.

Attention! Line 45 (`REF=~/Dropbox/phd/fastarefs/"$2"`), you will need to change the path to wherever your fasta reference is, for example: `REF=~/myrefs/"$2"`

About filtering; it keeps only reads

* whose length is > 140 bp (can change parameter line 4)
* whose Phred score is > 40 (can change parameter line 5)
* which have < 20% of their length soft-clipped (can change parameter line 76)

Script creates a folder **filterfastq** and puts it the fastq files converted back from the filtered bam.

To have a visual check of the aligment before vs after filtering, you can open them in IGV: load the right reference genome (Genomes > Load Genome From File...; pick your fasta sequence in folder *fastarefs*) and drag/drop bam file. There must be a bam.bai index file in the same folder and with the same name.

### 5– Run ampliCan

ampliCan is ran with the R scripts, eg. *amplican_slc24a5.R*
Before you run: create a *amplican* folder and you will need to change the paths.

ampliCan writes results files and reports (in *reports* folder) in the **amplican** folder.

For each sample, the two result files that are brought forward are: *config_summary.csv* and *events_filtered_shifted_normalized.csv* (in **alignments** folder).

## Analysis after ampliCan

All the *config_summary.csv* results files were collated with details about the loci and primers (from Supplementary file 1, sheets crRNAs_MiSeq and MiSeq_offtargets) into *MiSeq_amplicanresults.xlsx*.

Analysis was then done with different R scripts.

### Figure 2A & Figure 3B (right): frameshift stacked barplot

Script `MiSeq_frameshiftstack.R`, which takes *MiSeq_amplicanresults.xlsx* as input.

Also in `MiSeq_frameshiftstack.R` are computed different summary statistics reported in main text:
* Total number of reads in the dataset
* Total number of unique loci sequenced
* Total number of genes targeted + sequenced
* Total number of samples
* Coverage per sample: mean ± standard deviation
* % mutated reads: mean ± standard deviation
* % frameshifted reads: mean ± standard deviation

### Figure 1E,F and Figure 2D: proportion of frameshift alleles

Script `biallelicpropability.R`, which takes *MiSeq_amplicanresults.xlsx* as input.

### Figure 2F: off-targets

Script `MiSeq_offtarget.R`, which takes *MiSeq_amplicanresults.xlsx* as input.

### Figure 2B: indels similarity

Find the files in folder **events** (relates to the individual events, rather than the summary counts per sample).

For Figure 2B and C: all the *events_filtered_shifted_normalized.csv* results files were appended below each other as follows.

`appendCsv.command` is ran when in folder **allcsv**, which contains all the *events_filtered_shifted_normalized.csv* files. It creates file *events_all.csv*.

If multiple reads are exactly the same (including they contain the same mutation), ampliCan keeps track of it in the *counts* column in the *events_filtered_shifted_normalized.csv* files, as this avoids writing a new line to the file. However, the analysis here will use a different definition of 'duplicate', as any mutation that is the same. While it is probably okay to use *events_all.csv* directly, I felt it was safer to 'unfold' the *counts* first, so writing each individual mutation on its own line.

For example: if read A is exactly the same as read B and both have a 5-bp insertion, ampliCan will write the 5-bp insertion and keep track of the two reads by writing 2 in column *counts*. Instead, I will 'unfold' these counts: writing the 5-bp insertion from read A and the 5-bp insertion from read B on two separate lines. I then handle duplicates slightly differently.

Script `unfoldAllEvents.R` takes *events_all.csv* as input, 'unfolds' the *counts*, and writes file *events_all_unfolded.csv*.

For Figure 2B: script `howSimilarDifferentSamples.R`, which takes *events_all_unfolded.csv* as input.

### Figure 2C: histogram of indel lengths frequencies

Script `allevents_histo.R`, which takes *events_all.csv* as input.

### Main text: diversity in null alleles

Script `quantifyDiversity.R`, which takes *events_all_unfolded.csv* as input.

Also in `MiSeq_frameshiftstack.R` are computed different summary statistics reported in main text:
* Proportion of alleles with a frameshift mutation when targeting three loci.
* How much this proportion increased when targeting a fourth locus (for *slc24a5* and *tyr*)

---

Feel free to get in touch for questions

  * [![alt text][1.2]][1] [@francois_kroll](https://twitter.com/francois_kroll)

  * :email: francois@kroll.be

<!-- icons with padding -->
[1.1]: http://i.imgur.com/tXSoThF.png (twitter icon with padding)

<!-- icons without padding -->
[1.2]: http://i.imgur.com/wWzX9uB.png (twitter icon without padding)

<!-- links to your social media accounts -->
[1]: https://twitter.com/francois_kroll
