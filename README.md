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

  -B  Batchfile. If you have many sequence datasets to screen, a batchfile will save time. The first line of the batchfile contains specifications for the whole batch (-r -c -n -t options here with spaces in between arguments and values). From line 2 onwards, each line represents a sequence dataset with the first field specifying the basename, the second field specifying read 1 and the third field, if present, represents read 2 of the pair. Fields must be tab-separated. If no file is specified in field 3, then it is assumed to be a single end dataset.

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

## Batch mode

If you have many datasets to check, it may be convenient to run this analysis in batch mode.
Specify a batchfile using the `-B` argument:

```
./contam  -B batchfile.txt
apptainer run --bind $PWD:/contam contam.sif -B batchfile.txt 

```

The batch file has overall run options on the first line, including `-r -c -n -t` options separated by spaces.
From line 2 onwards the first field represents the basename of the dataset, field 2 is the first fastq file and
field 3 is optional second read in the pair.
Lines 2 onward need to be separated by tabs.
Here is an example.
 
```
-r ref -c contam_folder -n 1000000 -t 16
P1	P1_R1.fastq.gz	P1_R2.fastq.gz
P2	P2_R1.fastq.gz	P2_R2.fastq.gz
P2_single	P2_R1.fastq.gz
```

Batch runs generate the same output as single files, and for convenience, an additional file is created
(name:<batchfile>_result.tsv) to summarise all datasets in the batch. 
The values shown are ratio of contaminant:host unique reads expressed as a percentage.
Each column is a potential contaminant and each row is a sample.
