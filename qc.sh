# Basic quality control for mapping PE illumina data to a distant reference

# Rob Lanfear, December 2016

# A few things to set before you go

inputf="/disks/dacelo/data/raw_data/Project_SN7001117R_0083_CKulheim_LBronham_Melaleuca/"
outputbase="/disks/dacelo/data/QC/test/"
ref="/disks/dacelo/data/raw_data/active_refs/Emel.fa.gz" # reference file as a fasta
gff="/disks/dacelo/data/raw_data/active_refs/Egrandis_genes_chr1_to_chr11.gff3"
adaptors="/disks/dacelo/data/programs/bbmap/resources/adapters.fa"
threads=50 # number of threads to use
minlen=50 # minimum length of read to keep after trimming
trimq=0 # trim bases with quality < this 


# set up dirs
outputrawqc=$outputbase"rawqc/"
outputtrimqc=$outputbase"trimmedqc/"
outputtrimreads=$outputbase"trimmed_reads/"
ngmout=$outputbase"ngm/"    
indexcov=$outputbase"indexcov/"
mkdir $outputbase
mkdir $ngmout
mkdir $outputrawqc
mkdir $outputtrimqc
mkdir $outputtrimreads

# run bbduk on all pairs of samples
echo "Trimming with bbduk"

for in1 in $(find $inputf -name "*R1_001.fastq.gz"); do
    in2=${in1%%R1_001.fastq.gz}"R2_001.fastq.gz"
    echo "running bbduk on"
    echo $in1
    echo $in2

    f1=$(basename ${in1%%R1_001.fastq.gz}"R1_001_trimmed.fastq.gz")
    f2=$(basename ${in1%%R1_001.fastq.gz}"R2_001_trimmed.fastq.gz")

    out1=$outputtrimreads$f1
    out2=$outputtrimreads$f2

    sampleid=$outputtrimreads${f1%%R1_001_trimmed.fastq.gz}

    bbduk.sh in1=$in1 in2=$in2 out1=$out1 out2=$out2 minlen=$minlen k=25 mink=8 ktrim=r ref=$adaptors hdist=1 overwrite=f qtrim=rl trimq=$trimq t=$threads bhist=$sampleid"bhist.txt" qhist=$sampleid"qhist.txt" gchist=$sampleid"gchist.txt" aqhist=$sampleid"aqhist.txt" lhist=$sampleid"lhist.txt" > $sampleid"bbduk_log.txt"

done

# run fastqc on all the raw and trimmed data files
echo "Running fastqc"
find $inputf -name '*.fastq.gz' | xargs fastqc  -o $outputrawqc -t $threads
find $outputtrimreads -name '*.fastq.gz' | xargs fastqc  -o $outputtrimqc -t $threads

# map reads to E. grandis reference
# our trimmed files look like: RL41_S1_R1_001_trimmed.fastq.gz  RL41_S1_R2_001_trimmed.fastq.gz
echo "Mapping to reference with NGM"

for in1 in $(find $outputtrimreads -name "*R1_001_trimmed.fastq.gz"); do
    in2=${in1%%R1_001_trimmed.fastq.gz}"R2_001_trimmed.fastq.gz"
    echo "mapping files: "
    echo $in1
    echo $in2

    # output file setup
    f1=$(basename $in1)
    id=${f1%%R1_001_trimmed.fastq.gz}
    outsamngm=$ngmout$id"trimmed.sam"

    ngm -t $threads -p -r $ref -1 $in1 -2 $in2 -o $outsamngm
    outbamngm=$ngmout$id"trimmed.bam"
    samtools view -bS -@ $threads $outsamngm > $outbamngm
    samtools sort -@ $threads  $outbamngm -o $outbamngm
    rm $outsamngm
    
    echo "running qualimap on sorted bams"
    outqualimap_all=$ngmout$id"qualimap_all/"
    qualimap bamqc -bam $outbamngm -outdir $outqualimap_all -nt $threads -c
    outqualimap_gff=$ngmout$id"qualimap_gff/"
	qualimap bamqc -bam $outbamngm -outdir $outqualimap_gff -gff $gff -nt $threads -c

done

echo "indexing bams"
ls ${ngmout}*.bam | parallel "samtools index {}"

echo "running indexcov"
goleft indexcov --directory $indexcov --sex "" ${ngmout}"*.bam"

echo "running multiqc"
multiqc $outputbase -o $outputbase