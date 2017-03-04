# Quality control for short-read sequencing data from non-model organisms

## What?

This repository contains some scripts that can be adapted to perform basic quality control on short read (e.g. Illumina) data from non-model species. Specifically, we are interested in situations where one has to map paired-end short read data to a reference genome from another species. If your species has a reference genome, you should probably try Teaser first (https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0803-1)

## Why? 

The Lanfear lab does a lot of genome sequencing of Eucalyptus species. We have a great reference genome for Eucalyptus grandis, but we don't often work with data from E. grandis. So, we need to be careful to look in detail at the quality of our data and the performance of different mapping algorithms etc. 

## Getting started

Here, you can obviously substitute the reference genome for whatever the closest (or otherwise best) reference genome is to the species for which you have data. But from here on in I'll just assume that you, like me, work on Eucalyptus.

1. Go and get the latest version of the E. grandis genome assembly from here: https://phytozome.jgi.doe.gov/pz/portal.html#!info?alias=Org_Egrandis

2. Install the following software:

	* bbtools: http://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/installation-guide/
	* NextGenMap: http://cibiv.github.io/NextGenMap/
	* BWA: https://github.com/lh3/bwa
	* Stampy: http://www.well.ox.ac.uk/stampy
	* fastqc: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
	* multiqc: http://multiqc.info/
	* qualimap: http://qualimap.bioinfo.cipf.es/
	* samtools: https://github.com/lh3/samtools
	* goleft: https://github.com/brentp/goleft
	* GNU parallel: https://www.gnu.org/software/parallel/

3. Get a gff file (you could use a bed file too) of the region of the genome you are most interested in. For our analyses, we know that the Eucalyptus genomes have lots of repeats (~40%) and we know that we have no hope of using these regions for inference. So, as well as looking at how reads map to the whole genome, we also look at the subset of the genome that corresponds just to genes, here. This file comes with the E. grandis reference in the ```/annotation``` folder, and is called ```Egrandis_297_v2.0.gene.gff3.gz```. For my purposes, I am only interested in the first 11 scaffolds, which in the E. grandis genome correspond to 11 pretty well figured-out chromosomes. The other ~50k scaffolds are little pieces that were hard to assemble. I don't make inferences from those. Here's a clunky one-liner to get just the genes on the first 11 chromosomes (after unzipping the gff3.gz file):

```awk '$1=="Chr01" || $1=="Chr02" || $1=="Chr03" || $1=="Chr04" || $1=="Chr05" || $1=="Chr06" || $1=="Chr07" || $1=="Chr08" || $1=="Chr09" || $1=="Chr10" || $1=="Chr11"' genes.gff3 > genes_chr1_to_chr11.gff3```

4. Make sure that the appropriate programs are in your path by editing your ```./bashrc```. We had some issues with with GNU parallel at the start. A useful tip here is to just type ```parallel --version```. If you're not getting the GNU version (some linux distros have other defaults) edit your ```~/.basrc``` appropriately. 

## Overview

For the gory details, just look at the scripts. I'll provide a short overview here though. 

A disclaimer is that this code is neither particularly optimal nor particularly beautiful. But it's optimal enough and beautiful enough for us. I'm sharing it because I've benefitted so much from others sharing their own code that I hope this can be of some help (maybe only to my own students, but it's a start).

Also note that there are no commandline arguments. That's by design. Without commandline arguments, the script itself is a perfect record of what we did. This helps becuase when I get new data, I drop a copy of this script into the folder with the fastq files, then run that version of the script. Then, in a year when I can't remember exactly what I did for the QC, it's all there in the script. 

## compare.sh

### What it does

Compares a few common pieces of mapping software, by mapping a single sample of reads (you could easily subset your reads if mapping one sample takes too long). All we do is map a single sample to the reference with all the mappers, then use FastQC to look at how it went on the whole genome, and on the subset we're interested in (the genes, for us).

### Settings

* **R1**: the R1 reads file in fastq.gz format
* **R2**: the R2 reads file in fastq.gz format
* **outputbase**: the output file where you want all the processed data, QC reports, and summary stats to go. Output will be stored in subfolders of this
* **ref**: The reference genome. A quirk of my laziness is that I want you to put a copy of this in $outputbase yourself, and I'll work with it from there.
* **gff**: the gff (or bed) file that specifies a subset of the genome you're interested in

### Running it

I just run it like this:

```sh compare.sh &> log.txt```

You can extract the time each analysis took from the frequent date stamps. I just do this by hand - it's only 4 mappers. 

### Output

A folder for each aligner (e.g. ```/ngm```), which contains:

1. The .bam and .bai files
2. A copy of the reference genome, and whatever index the mapper created for it
3. ```/qualimap_all``` qualimap results from looking at the whole .bam file
4. ```/qualimap_gff``` qualimap results from looking at the subset of the .bam file specified by the input .gff file

### What to look at

Compare the qualimap results for all of the aligners. Beware that the mapping qualities are not comparable - each aligner calculates mapping qualities in a different (sometimes very different) way. Particular things to pay attention to are

1. Is the error rate roughly what you expect (i.e. the % divergence of your species from the reference you're using)
2. Is the coverage on the mappable regions (the bits you're interested in in the gff file) roughly what you expected given your sequencing setup (i.e. given the GB of data you have and the genome size of your species, but be careful if the ref. you're using has a very different genome size from your focal spp)
3. Do the insert sizes correspond to what you expect from the fragment sizes you input (e.g. compare to your bioanalyser results from pre-sequencing)?
4. Which mapper maps the most reads? 
5. How's the coverage across the chromosomes? 
6. How's the insert size across the chromosomes?

## qc.sh

### What it does

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

Why am I using bbduk? I don't know. It would be better to just use whatever the best mapper you identified from ```compare.sh``` I just haven't got round to updating the script yet. For me the best mapper for this is clearly NextGenMap - it's very accurate and very fast. 

### Settings

* **inputf**: the input file where all your reads are (they can be in subfolders)
* **outputbase**: the output file where you want all the processed data, QC reports, and summary stats to go. Output will be stored in subfolders of this
* **bbduk**: the location of the bbduk.sh script (from when you installed bbtools)
* **bbmap**: the location of the bbmap.sh script (from when you installed bbtools)
* **qualimap**: the location of the qualimap executable
* **adaptors**: the adaptor file you want to use with bbduk (we usually just use the adapters.fa file, which contains most commonly used adaptors.
* **threads**: maximum number of threads to use
* **minlen**: minimum length of read to keep after quality and adaptor trimming
* **trimq**: trim bases with quality < this 

### Output

* ```/bbmap``` sorted BAMs and various summary outputs from mapping with bbmap
* ```/rawqc``` fastqc results on the raw sequencing reads
* ```/trimmedqc``` fastqc results on the trimmed sequencing reads
* ```/trimmed_reads``` trimming stats from bbduk (lots of text files)
* ```/multiqc_data``` data for the multiQC report
* ```/indexcov``` information from indexcov
* ```multiqc_report.html``` the multiQC report itself

### What to look at

1. Use the multiQC report for a very useful (but necessarily rather limited) overall summary of the data. This is the first thing we look at to detect issues. 
2. Examine the stats from bbduk in ```/trimmed_reads```. This will tell you useful information about adaptors and quality trimming.
3. Examine the fastqc.html files in ```/rawqc``` and ```/trimmedqc```
4. Look at the indexcov HTML report for any issues, particularly for differences between samples that might indicate problems with extraction / sequencing. 
