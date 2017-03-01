# Quality control for short-read sequencing data from Eucalyptus species

## What?

This repository contains some scripts that can be adapted to perform basic quality control on short read (e.g. Illumina) data from Eucalyptus species.

## Why? 

Because we do a lot of that kind of sequencing, and it's a pain figuring this stuff out from scratch. If you have suggestions, we really would love to hear them. Please leave suggestions in the issues. 

## Getting started

1. Go and get the latest version of the E. grandis genome assembly from here: https://phytozome.jgi.doe.gov/pz/portal.html#!info?alias=Org_Egrandis

2. Install the following software:

	* bbtools: http://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/installation-guide/
	* fastqc: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
	* multiqc: http://multiqc.info/
	* qualimap: http://qualimap.bioinfo.cipf.es/
	* samtools: https://github.com/lh3/samtools
	* goleft: https://github.com/brentp/goleft
	* GNU parallel: https://www.gnu.org/software/parallel/

Make sure that GNU parallel is in your path such that typing ```parallel``` calls it. If you want to know which version you call, just type ```parallel --version```. If you're not getting the GNU version (some linux distros have other defaults) edit your ```~/.basrc``` appropriately. 

## Overview

For the gory details, just look at the scripts. I'll provide a short overview here though. 

A disclaimer is that this code is neither particularly optimal nor particularly beautiful. But it's optimal enough and beautiful enough for us. I'm sharing it because I've benefitted so much from others sharing their own code that I hope this can be of some help (maybe only to my own students, but it's a start).

Also note that there are no commandline arguments. That's by design. Without commandline arguments, the script itself is a perfect record of what we did. This helps becuase when I get new data, I drop a copy of this script into the folder with the fastq files, then run that version of the script. Then, in a year when I can't remember exactly what I did for the QC, it's all there in the script. 

### qc.sh

#### What it does

Performs the basic QC on the reads, in the following steps:

1. Uses *bbduk* for adaptor then quality trimming of each set of reads
2. Uses *fastqc* to make some basic measurements of the raw and trimmed reads
3. Uses *bbmap* to map the trimmed reads to the E. grandis genome
4. Uses *samtools* to turn the SAM files to sorted BAM files
5. Deletes all the large files (the SAMs and the trimmed read files)*
6. Uses *samtools* and GNU parallel to index all the bams
7. Uses *indexcov* to visualise various coverage stats across all samples
6. Uses *multiQC* to summarise most of the QC data

*It would be easy not to delete these files of course, but since the only point of this script is to run QC, we usually don't want to keep the trimmed reads. We can always re-trim the files for the final analysis if we like the settings we used.

#### Settings

* **inputf**: the input file where all your reads are (they can be in subfolders)
* **outputbase**: the output file where you want all the processed data, QC reports, and summary stats to go. Output will be stored in subfolders of this
* **bbduk**: the location of the bbduk.sh script (from when you installed bbtools)
* **bbmap**: the location of the bbmap.sh script (from when you installed bbtools)
* **qualimap**: the location of the qualimap executable
* **adaptors**: the adaptor file you want to use with bbduk (we usually just use the adapters.fa file, which contains most commonly used adaptors.
* **threads**: maximum number of threads to use
* **minlen**: minimum length of read to keep after quality and adaptor trimming
* **trimq**: trim bases with quality < this 

#### Output

* ```/bbmap``` sorted BAMs and various summary outputs from mapping with bbmap
* ```/rawqc``` fastqc results on the raw sequencing reads
* ```/trimmedqc``` fastqc results on the trimmed sequencing reads
* ```/trimmed_reads``` trimming stats from bbduk (lots of text files)
* ```/multiqc_data``` data for the multiQC report
* ```/indexcov``` information from indexcov
* ```multiqc_report.html``` the multiQC report itself

#### What to look at

1. Use the multiQC report for a very useful (but necessarily rather limited) overall summary of the data. This is the first thing we look at to detect issues. 
2. Examine the stats from bbduk in ```/trimmed_reads```. This will tell you useful information about adaptors and quality trimming.
3. Examine the fastqc.html files in ```/rawqc``` and ```/trimmedqc```
4. Look at the indexcov HTML report for any issues, particularly for differences between samples that might indicate problems with extraction / sequencing. 
