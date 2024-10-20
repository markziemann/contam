# contam

`contam` is a bash script that uses BWA MEM alignment to quantify contamination from
Illumina sequences.
This is particularly useful for screening RNA-seq data derived from in vitro cultures where mycoplasma
contamination is fairly common.
It is tested to work with single- and paired-end RNA-seq and DNA-seq, and might work for other sequence data types.
An apptainer image is provided to allow portability.

## Installation - native

Download the 'contam' shell script.
Move it to the working directory and enable execution with `chmod +x contam`.

Dependancies: 'skewer', 'bwa' and 'samtools', which are available through `apt` on Debian-based
systems.

execution: `./contam -h`

## Installation - Apptainer

Apptainer is a containerisation framework similar to Docker, but it does not require root permissions
for operation.
Apptainer pre built packages are available; see the [docs](https://apptainer.org/docs/admin/main/installation.html#install-from-pre-built-packages) for details.

After apptainer is installed, you can download and run the apptainer image file (`contam.sif`).

```
wget https://ziemann-lab.net/public/contam/contam.sif
apptainer run --bind $PWD:/contam contam.sif -h
```

## Usage and options

From the help page:

```
contam is a script for the quantification of contaminant sequences in Illumina short read sequence data including genomic and transcriptome samples.

Usage: contam [-h] [-v] <-r HOST_REFERENCE_GENOME_FOLDER> <-c CONTAMINANTS_FOLDER> [-n NUM_READS] [-t NUM_THREADS] <-1 FASTQ_FILE_1> <-2 FASTQ_FILE_2> <-b FASTQ_BASE_NAME>

  -h  Help. Display this message and quit.

  -v  Version. Print version number and quit.

  -r  Reference sequence (host) folder. Must contain one gzip compressed fasta file with '.fa.gz' suffix.

  -c  Contaminant sequence folder. May contain one or many compressed fasta files.

  -n  Number of read pairs to analyze. Default is 1000000.

  -t  Number of parallel threads. Default is 2.

  -1  First fastq file corresponding to read 1 of the pair.

  -2  Second fastq file corresponding to read 2 of the pair. Omit for single-end sequencing runs.

  -b  Basename of the sequence dataset. This is the name of the folder where results will be saved.

Output: For each sequence dataset, a folder is created which contains the intermediate files
and final results, called "result.tsv".
The result.tsv file contains four columns, the contaminant name (abbreviated to 15 characters),
The number of reads that mapped to the host genome only, the number of reads that mapped to
the contaminant genome only and the number of reads that mapped to both.
Each line in the file corresponds to a different contaminant.

```

### Native usage

```
./contam -r ref -c contam_folder -n 10000 -t 8 -1 test_R1.fq.gz -2 test_R2.fq.gz -b test
```

## Apptainer usage

```
apptainer run --bind $PWD:/contam contam.sif -r ref -c contam_folder -n 10000 -t 8 -1 test_R1.fq.gz -2 test_R2.fq.gz -b test
```
