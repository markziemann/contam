#!/bin/bash

REF=$1
CONTAM=$2
NREADS=$6
THREADS=$7
FQZ1=$4
FQZ2=$5
BASE=$3
FQ1=$BASE/R1.fastq
FQ2=$BASE/R2.fastq
FQT1=$BASE/R1-trimmed-pair1.fastq
FQT2=$BASE/R1-trimmed-pair2.fastq

## CHECKS
if [ ! -d "$REF" ] ; then
  echo "'ref' directory does not exist."
  exit
fi

if [ ! -d "$CONTAM" ] ; then
  echo "'contam' directory does not exist."
  exit
fi

if [ ! -r "$FQZ1" ] ; then
  echo "'FQ1' does not exist."
  exit
fi

if [ ! -r "$FQZ2" ] ; then
  echo "'FQ2' does not exist."
  exit
fi

SFX1=$(echo $FQZ1 | egrep -c  '(fastq.gz$|fa.gz$)' )
if [ $SFX1 -eq 0 ] ; then
  echo FQ1 needs to have a fa.gz or fastq.gz suffix
fi

SFX2=$(echo $FQZ2 | egrep -c  '(fastq.gz$|fa.gz$)' )
if [ $SFX2 -eq 0 ] ; then
  echo FQ2 needs to have a fa.gz or fastq.gz suffix
fi

# set ref seq
REFCNT=$(find "$REF" | grep -c fa.gz$)
if [ $REFCNT -ne 1 ] ; then
  echo "There needs to be one reference sequence with a 'fa.gz' suffix!"
  exit
fi

REFFA=$(find "$REF" | grep fa.gz$ | head -1)

BWT=$REFFA.bwt

if [ ! -r $BWT ] ; then
  bwa index $REFFA
fi

# set up contam sequences
CONTAMCNT=$(find "$CONTAM" | grep -c fa.gz$)
if [ $CONTAMCNT -eq 0 ] ; then
  echo "There needs to be at least one contaminant genome sequence with a 'fa.gz' suffix!"
  exit
fi

for CONTAMFA in $CONTAM/*fa.gz ; do
  bwa index $CONTAMFA
done

## start the analysis
mkdir $BASE
NLINES=$(echo $NREADS | awk '{print $1 * 4}')

zcat $FQZ1 | head -$NLINES > $FQ1
zcat $FQZ2 | head -$NLINES > $FQ2

skewer -t $THREADS -q 20 $FQ1 $FQ2 2> /dev/null

RES=$BASE/result.tsv

HOSTCNT=$(bwa mem -t $THREADS $REFFA $FQT1 $FQT2 2> /dev/null \
| samtools view -q 10 -F 4 - 2> /dev/null \
| awk '!arr[$1]++' | wc -l )

for CONTAMFA in $CONTAM/*fa.gz ; do

  >&2 echo
  >&2 echo "Starting analysis of $FQZ1 and $FQZ2 with contaminant $CONTAMFA"
  >&2 echo

  SHORT=$(echo $CONTAMFA | rev | cut -f1 -d '/' | rev | cut -c-15 | sed 's/.fa.gz//')
  CONTAMBAM=$BASE/$SHORT.bam
  MAPPEDFQ1=$CONTAMBAM.R1.fastq
  MAPPEDFQ2=$CONTAMBAM.R2.fastq
  MAPPEDFQSINGLES=$CONTAMBAM.SINGLES.fastq

  bwa mem -t $THREADS $CONTAMFA $FQT1 $FQT2 2> /dev/null \
  | samtools view -q 10 -F 4 -uSh - 2> /dev/null \
  | samtools sort -o $CONTAMBAM - 2> /dev/null

  samtools index $CONTAMBAM

  CONTAMCNT=$(samtools view -q 10 -F 4 $CONTAMBAM | awk '!arr[$1]++' | wc -l )

  if [ $CONTAMCNT -gt 0 ] ; then
    samtools fastq -F 4 -1 $MAPPEDFQ1 -2 $MAPPEDFQ2 -s $MAPPEDFQSINGLES $CONTAMBAM 2> /dev/null

    REFMAPCNT=$(bwa mem -t $THREADS $REFFA $MAPPEDFQ1 $MAPPEDFQ2 2> /dev/null \
    | samtools view -q 10 -F 4 - 2> /dev/null \
    | awk '!arr[$1]++' | wc -l )

  else
    REFMAPCNT=0
  fi

  echo $SHORT $HOSTCNT $CONTAMCNT $REFMAPCNT

  >&2 echo
  >&2 echo "Done!"
  >&2 echo

done | tr ' ' '\t' >> $RES

awk '{OFS="\t"} {print $1,$2-4,$3-$4,$4}' $RES > $RES.tmp

echo CONTAMINANT HOST_ONLY CONTAM_ONLY BOTH | tr ' ' '\t' > $RES
cat $RES.tmp >> $RES
rm $RES.tmp

cat $RES

