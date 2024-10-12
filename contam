#!/bin/bash
VERSION=0.1

# Parse arguments
usage() {
    echo
    echo "$0 is a script for the quantification of contaminant sequences in Illumina short read sequence data including genomic and transcriptome samples."
    echo
    echo "Usage: $0 [-h] [-v] <-r HOST_REFERENCE_GENOME_FOLDER> <-c CONTAMINANTS_FOLDER> [-n NUM_READS] [-t NUM_THREADS] <-1 FASTQ_FILE_1> <-2 FASTQ_FILE_2> <-b FASTQ_BASE_NAME>"
    echo
    echo "  -h  Help. Display this message and quit."
    echo "  -v  Version. Print version number and quit."
    echo "  -r  Reference sequence (host) folder. Must be one gzip compressed fasta file with '.fa.gz' suffix."
    echo "  -c  Contaminant sequence folder. May contain one or many compressed fasta files."
    echo "  -n  Number of read pairs to analyze. Default is 1000000."
    echo "  -t  Number of parallel threads. Default is 2."
    echo "  -1  First fastq file corresponding to read 1 of the pair."
    echo "  -2  Second fastq file corresponding to read 2 of the pair."
    echo "  -b  Basename of the sequence dataset. This is the name of the folder where results will be saved."
    echo
    exit
}

if [ $# -eq 0 ] ; then
  usage
fi

NREADS=100000
THREADS=2

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -h|--help) usage ; exit ;;
        -v|--version) echo "Version $VERSION" ; exit ;;
        -r|--ref) REF="$2" ; shift ;;
        -c|--contaminants) CONTAM="$2" ; shift ;;
        -1|--fastq1) FQZ1="$2" ; shift ;;
        -2|--fastq2) FQZ2="$2" ; shift ;;
        -n|--numreads) NREADS="$2" ; shift ;;
        -t|--threads) THREADS="$2" ; shift ;;
        -b|--basename) BASE="$2" ; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

FQ1=$BASE/R1.fastq
FQ2=$BASE/R2.fastq
FQT1=$BASE/R1-trimmed-pair1.fastq
FQT2=$BASE/R1-trimmed-pair2.fastq

## CHECKS
if [ -z "$REF" ] ; then
  echo "Error: HOST_REFERENCE_GENOME_FOLDER not set."
  exit
fi

if [ -z "$CONTAM" ] ; then
  echo "Error: CONTAMINANTS_FOLDER not set."
  exit
fi

if [ -z "$FQZ1" ] ; then
  echo "Error: 'FASTQ_FILE_1' not set."
  exit
fi

if [ -z "$FQZ2" ] ; then
  echo "Error: 'FASTQ_FILE_2' not set."
  exit
fi

if [ ! -d "$REF" ] ; then
  echo "HOST_REFERENCE_GENOME_FOLDER does not exist."
  exit
fi

if [ ! -d "$CONTAM" ] ; then
  echo "Error: 'CONTAMINANTS_FOLDER' does not exist."
  exit
fi

if [ ! -r "$FQZ1" ] ; then
  echo "Error: 'FASTQ_FILE_1' does not exist."
  exit
fi

if [ ! -r "$FQZ2" ] ; then
  echo "Error: 'FASTQ_FILE_2' does not exist."
  exit
fi

if [ -z "$BASE" ] ; then
  echo "Error: FASTQ_BASE_NAME not set."
  exit
fi

if [ "$THREADS" -gt 32 ] ; then
  echo "The number of parallel threads is capped at 32."
  THREADS=32
fi

SFX1=$(echo $FQZ1 | egrep -c  '(fastq.gz$|fa.gz$)' )
if [ $SFX1 -eq 0 ] ; then
  echo FASTQ_FILE_1 needs to have a fa.gz or fastq.gz suffix
fi

SFX2=$(echo $FQZ2 | egrep -c  '(fastq.gz$|fa.gz$)' )
if [ $SFX2 -eq 0 ] ; then
  echo FASTQ_FILE_2 needs to have a fa.gz or fastq.gz suffix
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
  echo "Indexing the reference genome. This might take a while"
  bwa index $REFFA 2> /dev/null
else
  echo "Reference genome is already indexed."
fi

# set up contam sequences
CONTAMCNT=$(find "$CONTAM" | grep -c fa.gz$)
if [ $CONTAMCNT -eq 0 ] ; then
  echo "There needs to be at least one contaminant genome sequence with a 'fa.gz' suffix!"
  exit
fi

for CONTAMFA in $CONTAM/*fa.gz ; do
  if [ ! -r $CONTAMFA.bwt ] ; then
    echo "Indexing contaminant genome $CONTAMFA"
    bwa index $CONTAMFA 2> /dev/null
  else
    echo "$CONTAMFA already indexed."
  fi
done

## start the analysis
mkdir $BASE 2> /dev/null

NLINES=$(echo $NREADS | awk '{print $1 * 4}')

echo "Starting fastq processing now."

zcat $FQZ1 | head -$NLINES > $FQ1
zcat $FQZ2 | head -$NLINES > $FQ2

skewer -t $THREADS -q 20 $FQ1 $FQ2 2> /dev/null > /dev/null

RES=$BASE/result.tsv
>$RES

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

cat $RES | column -t

echo
