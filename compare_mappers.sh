# Compare BWA-MEM, bbmap, nextgenmap, stampy for mapping to a distant reference

# a little shell script to compare a few mappers. 
# note: I drop the reference.fa.gz file into the base directory ($outputbase) below
# also note: I do some fairly ill-advised things like cd-ing around directories

# Rob Lanfear, March 2016

# forward and reverse input reads
R1="/disks/dacelo/data/raw_data/Project_SN7001117R_0083_CKulheim_LBronham_Melaleuca/Sample_M1a_index2/M1a_index2_CGATGTAT_L003_R1_001.fastq.gz"
R2="/disks/dacelo/data/raw_data/Project_SN7001117R_0083_CKulheim_LBronham_Melaleuca/Sample_M1a_index2/M1a_index2_CGATGTAT_L003_R2_001.fastq.gz"

# directory where you want to do stuff
outputbase="/disks/dacelo/data/compare_mappers/"

# reference genome
ref="Emel.fa.gz" # reference file as a fasta

# gff file defining subset of reference genome - qualimap will run on both the whole genome and on this subset
# sensible to define a subset that you might use for your final analyses
gff=$outputbase"Egrandis_genes_chr1_to_chr11.gff3"

threads=20 # number of threads to use

# set up dirs
outBWAMEM=$outputbase"BWAMEM/"
outbbmap=$outputbase"bbmap/"
outngm=$outputbase"ngm/"
outstampy=$outputbase"stampy/"
mkdir $outputbase
mkdir $outBWAMEM
mkdir $outbbmap
mkdir $outngm
mkdir $outstampy


cd $outngm
cp ../$ref $ref
date
echo "indexing and mapping with ngm"
time ngm -t $threads -p -r $ref -1 $R1 -2 $R2 -o out.sam 
date
echo "indexing and mapping with ngm done"
echo "converting sams to bams"
time samtools view -bS -@ $threads out.sam > out.bam
date
echo "sorting and indexing bams"
time samtools sort -@ $threads out.bam -o out.bam
time samtools index out.bam
date
echo "running qualimap"
time qualimap bamqc -bam out.bam -outdir $outngm"qualimap_all/" -nt $threads -c
time qualimap bamqc -bam out.bam -outdir $outngm"qualimap_gff/" -gff $gff -nt $threads -c
date
echo "Done mapping with ngm"


cd $outbbmap
cp ../$ref $ref
date
echo "indexing and mapping with bbmap"
time bbmap.sh in1=$R1 in2=$R2 ref=$ref out=out.sam path=$outbbmap t=$threads
date
echo "indexing and mapping with bbmap done"
echo "converting sams to bams"
time samtools view -bS -@ $threads out.sam > out.bam
date
echo "sorting and indexing bams"
time samtools sort -@ $threads out.bam -o out.bam
time samtools index out.bam
date
echo "running qualimap"
time qualimap bamqc -bam out.bam -outdir $outbbmap"qualimap_all/" -nt $threads -c
time qualimap bamqc -bam out.bam -outdir $outbbmap"qualimap_gff/" -gff $gff -nt $threads -c
date
echo "Done mapping with bbmap"

cd $outBWAMEM
cp ../$ref $ref
date
echo "indexing and mapping with BWA MEM"
echo "Creating index with BWA"
time bwa index $ref
date
echo "Mapping with BWA MEM"
time bwa mem $ref $R1 $R2 -t $threads > out.sam
date
echo "indexing and mapping with BWA MEM done"
echo "converting sams to bams"
time samtools view -bS -@ $threads out.sam > out.bam
date
echo "sorting and indexing bams"
time samtools sort -@ $threads out.bam -o out.bam
time samtools index out.bam
date
echo "running qualimap"
time qualimap bamqc -bam out.bam -outdir $outBWAMEM"qualimap_all/" -nt $threads -c
time qualimap bamqc -bam out.bam -outdir $outBWAMEM"qualimap_gff/" -gff $gff -nt $threads -c
date
echo "Done mapping with BWA MEM"

cd $outstampy
cp ../$ref $ref
date
echo "indexing and mapping with Stampy"
echo "indexing with Stampy"
time stampy.py -G index $ref
time stampy.py -g index -H index
date
echo "mapping with Stampy"
time stampy.py -g index -h index -t $threads -M $R1 $R2 > out.sam
date
echo "indexing and mapping with Stampy done"
echo "converting sams to bams"
time samtools view -bS -@ $threads out.sam > out.bam
date
echo "sorting and indexing bams"
time samtools sort -@ $threads out.bam -o out.bam
time samtools index out.bam
date
echo "running qualimap"
time qualimap bamqc -bam out.bam -outdir $outstampy"qualimap_all/" -nt $threads -c
time qualimap bamqc -bam out.bam -outdir $outstampy"qualimap_gff/" -gff $gff -nt $threads -c
date
echo "Done mapping with stampy"
