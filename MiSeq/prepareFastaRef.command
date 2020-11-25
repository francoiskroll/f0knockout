#!/bin/bash

shopt -s nullglob

echo
echo
echo
echo "---- [ CREATING REFERENCE FOR "$" $1 ] ----"
echo
echo
echo

REF=~/Dropbox/phd/fastarefs/"$1"
echo "$REF"

bwa index -a bwtsw "$REF"

samtools faidx "$REF"

PRE="$(echo "$1" | cut -d'.' -f 1)"
echo "$PRE"

DICT="$(echo $PRE$".dict")"
echo "$DICT"

DICTPATH=~/Dropbox/phd/fastarefs/"$DICT"
echo "$DICTPATH"

picard CreateSequenceDictionary \
  REFERENCE="$REF" \
  OUTPUT="$DICTPATH"

shopt -u nullglob
