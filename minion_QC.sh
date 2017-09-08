# Basic quality control for mapping MinION data to a reference

# Rob Lanfear and Miriam Schalamun

# A few things to set before you go
inputf="/data/nanopore/nanopore/RB7_A2/"
outputbase=$inputf
ref="/data/active_refs/Epau.fa.gz" # reference file as a fasta
gff="/data/active_refs/Egrandis_genes_chr1_to_chr11.gff3"
threads=25 # number of threads to use
mem_size='50G' # memory size for Qualimap
flowcellID="FLO-MIN107"
kitID="SQK-LSK108"

mkdir $outputbase

# basecall with albacore 
# the --disable-filtering is new for albacore 2.x, and means that we keep all the reads in one folder
# this allows us to look at all the data more easily, and determine our own cutoff flowcell-by-flowcell
# -q 0 puts all the reads into a single fastq file
read_fast5_basecaller.py -i $inputf -t $threads -s $outputbase -f $flowcellID -k $kitID -r -o fastq -q 0 --disable_filtering

# cat together the fastq's to one big fastq
fastq_file=$outputbase"reads.fastq"
cat $outputbase/workspace/*.fastq > $fastq_file

# minion QC
echo "Running MinION QC R script"
outmqc=$outputbase/"minionQC"
seqsum=$outputbase/"sequencing_summary.txt"
minion_QC.R $seqsum $outmqc

# map with ngmlr
outngmlr=$outputbase"ngmlr/"
mkdir $outngmlr
cd $outngmlr
echo "Mapping with ngmlr"
date
time ngmlr -t $threads -r $ref -q $fastq_file -o out.sam -x ont
echo "Done Mapping with ngmlr"
date
samtools view -bS -@ $threads out.sam > out.bam
samtools sort -@ $threads out.bam -o out.bam
samtools index out.bam
rm out.sam
qualimap bamqc -bam out.bam -outdir $outputbase"qualimap_all/" -nt $threads -c --java-mem-size=$mem_size
qualimap bamqc -bam out.bam -outdir $outputbase"qualimap_gff/" -gff $gff -nt $threads -c --java-mem-size=$mem_size

# stats on reads > various length (thanks to @gringer here: https://bioinformatics.stackexchange.com/questions/678/get-the-mapping-statistics-of-a-single-read-from-a-bam-file/696#696)
outbam=$outngmlr"out.bam"
samtools view -h $outbam |     awk -F'\t' '{if((/^@/) || (length($10)>1000)){print $0}}' |  samtools stats   | grep '^SN' | cut -f 2- > stats_1k.txt
samtools view -h $outbam |     awk -F'\t' '{if((/^@/) || (length($10)>2000)){print $0}}' |  samtools stats   | grep '^SN' | cut -f 2- > stats_2k.txt
samtools view -h $outbam |     awk -F'\t' '{if((/^@/) || (length($10)>10000)){print $0}}' | samtools stats   | grep '^SN' | cut -f 2- > stats_10k.txt
samtools view -h $outbam |     awk -F'\t' '{if((/^@/) || (length($10)>20000)){print $0}}' | samtools stats   | grep '^SN' | cut -f 2- > stats_20k.txt
samtools view -h $outbam |     awk -F'\t' '{if((/^@/) || (length($10)>100000)){print $0}}' | samtools stats  | grep '^SN' | cut -f 2- > stats_100k.txt
samtools view -h $outbam |     awk -F'\t' '{if((/^@/) || (length($10)>200000)){print $0}}' | samtools stats  | grep '^SN' | cut -f 2- > stats_200k.txt 

echo "Assembling with minimap & miniasm"
miniasmout=$outputbase"miniasm/"
mkdir $miniasmout
paf=$miniasmout"reads.paf.gz"
minimap -Sw5 -L100 -m0 -t 35 $fastq_file $fastq_file | gzip -1 > $paf
miniasm -f $fastq_file $paf > $miniasmout"miniasm.gfa"
awk '$1=="S" {print ">"$2"\n"$3} ' $miniasmout"miniasm.gfa" > $miniasmout"miniasm.fa"

echo "Assessing assembly with quast"
quastout=$outputbase"quast/"
mkdir $quastout
python quast.py  -t $threads -o $quastout --gene-finding --eukaryote -R $ref -G $gff $miniasmout"miniasm.fa"
