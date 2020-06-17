#!/bin/bash

# filtering parameters
declare -i length=140
declare -i phred=40
# maximum soft-clipping is 20% (below, cannot find how to pass the variable)

shopt -s nullglob

################################################################################

echo
echo
echo "---- [ ALIGNING"$" $1 TO"$" $2 ] ----"
echo
echo

FILE="$1"

# CHECK I HAVE FWD READS. IF OK: GET NAME OF RVS READS
DIR="$(echo "$1" | cut -d'_' -f 4)" # DIR should be R1 = Forward reads or R2 = Reverse reads

if [ "$DIR" = "R1" ]; then # do I have the Forward reads
# if yes: need to find the Reverse reads
  FWD="$FILE"

  HALF1="$(echo "$FWD" | cut -d'R' -f 1)" # everything before 'R'
  HALF2="$(echo "$FWD" | cut -d'_' -f 5)" # everything after '_'
  RVS="$(echo $HALF1$"R2_"$HALF2)" # based on the name of the Forward file, this should be the name of the Reverse file
else
  echo "Error: give me Forward reads"
fi

WELL="$(echo "$FWD" | cut -d'_' -f 1 | cut -d'-' -f 2)" # first cut gets everything before first '_'; second cut gets everything after first '-'

echo
echo "---- [ Found reads from well"$" $WELL ] ----"
echo

# alignment & convert/sort sam to bam
BAM="$(echo "$WELL"$".bam")"

# reference should be indexed with prepareFastaRef.command; all in ~/Dropbox/phd/fastarefs

REF=~/Dropbox/phd/fastarefs/"$2"

echo
echo "---- [ Searching reference at"$" $REF ] ----"
echo

echo
echo "---- [ Aligning"$" $FWD and"$" $RVS to"$" $REF ] ----"
echo

bwa mem -M -t 16 "$REF" "$FWD" "$RVS" \
  | samtools sort > "$BAM"

# index the sorted bam
samtools index "$BAM"

################################################################################

echo
echo
echo "---- [ FILTERING"$" $BAM ] ----"
echo
echo

# remove any reads shorter than length
FILTER1="$(echo "$WELL"$"_1.bam")"
samtools view -h "$BAM" | awk -v lgth="$length" 'length($10) > lgth || $1 ~ /^@/' | awk -v phr="$phred" '$5 > phr || $1 ~ /^@/' | samtools view -bS - > "$FILTER1"

# remove any reads with > 20% soft-clipped
FILTER2="$(echo "$WELL"$"_filter.bam")"

java -jar ~/packages/jvarkit/dist/samjdk.jar -e 'return record.getReadUnmappedFlag() || record.getCigar()==null || record.getCigar().getCigarElements().stream().filter(C->C.getOperator().isClipping()).mapToInt(C->C.getLength()).sum() / (double)record.getCigar().getReadLength() < 0.20;' "$FILTER1" \
  | samtools view -S -b > "$FILTER2"
samtools index "$FILTER2" # so final filtered bam can be opened in IGV

rm "$FILTER1"

################################################################################

echo
echo
echo "---- [ CONVERTING BACK TO FASTQ"$" $FILTER2 ] ----"
echo "---- [ INTO"$" $FWD ] ----"
echo
echo

BAMFILTER="$FILTER2"
FWDOUT="$FWD"

FILTER3="$(echo "$WELL"$"_3.bam")"

# generate right name of fastq files
FWDNAME="$(echo "$FWDOUT" | cut -d'.' -f 1)" # name of the file without .fastq.gz
FWDFQ="$(echo "$FWDNAME"$".fastq")"

HALF1="$(echo "$FWDOUT" | cut -d'R' -f 1)" # everything before 'R'
HALF2="$(echo "$FWDOUT" | cut -d'_' -f 5)" # everything after '_'
RVS="$(echo $HALF1$"R2_"$HALF2)" # based on the name of the Forward file, this should be the name of the Reverse file

RVSNAME="$(echo "$RVS" | cut -d'.' -f 1)"
RVSFQ="$(echo "$RVSNAME"$".fastq")"

samtools sort -n "$BAMFILTER" -o "$FILTER3"

mkdir filterfastq # fastq will have exact same names as original so put them in a separate folder
bedtools bamtofastq -i "$FILTER3" -fq  $"filterfastq/""$FWDFQ" -fq2 $"filterfastq/""$RVSFQ"
gzip -f $"filterfastq/""$FWDFQ"
gzip -f $"filterfastq/""$RVSFQ"

# delete unnecessary files
rm "$FILTER3"

################################################################################

shopt -u nullglob
