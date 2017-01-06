# Basic quality control for Eucalyptus illumina data

# Rob Lanfear, December 2016

# A few things to set before you go

inputf="/disks/dacelo/data/raw_data/Project_SN7001117R_0083_CKulheim_LBronham_Melaleuca/"
outputbase="/disks/dacelo/data/QC/Project_SN7001117R_0083_CKulheim_LBronham_Melaleuca/"
ref="/disks/dacelo/data/raw_data/indices/Emel.fa.gz" # reference file as a fasta
bbduk="/disks/dacelo/data/programs/bbmap/bbduk.sh"
bbmap="/disks/dacelo/data/programs/bbmap/bbmap.sh"
goleft="/disks/dacelo/data/programs/goleft/goleft_linux64"
qualimap="/disks/dacelo/data/programs/qualimap_v2.2.1/qualimap"
adaptors="/disks/dacelo/data/programs/bbmap/resources/adapters.fa"

# for more control of the parameters, you can change them in the commandlines below
# of course, if there's a parameter you really want to explore, make it a variable here.

threads=50 # number of threads to use
minlen=50 # minimum length of read to keep after trimming
trimq=20 # trim bases with quality < this 

# chromosomes on which to run indexcov
# the default is the first 11 (which are the good nuclear chromosomes of E. grandis)
# and the last one (which by our own convention is the choloroplast genome)
nuc_chroms=$(zcat $ref | grep -m11 '>' | cut -c 2- )
chl_chrom=$(zcat $ref | tac | grep -m1 '>' | cut -c 2- )
chroms = $nuc_chroms' '$chl_chrom


echo "indexcov will run on the following chromosomes (not now, later...)"
echo $chroms

outputrawqc=$outputbase"rawqc/"
outputtrimqc=$outputbase"trimmedqc/"
outputtrimreads=$outputbase"trimmed_reads/"
indexcov=$outputbase"indexcov/"

echo $outputbase
echo $outputrawqc
echo $outputtrimqc
echo $outputtrimreads
echo $indexcov

# set up dirs
mkdir $outputbase
mkdir $outputbase"bbmap"
mkdir $outputrawqc
mkdir $outputtrimqc
mkdir $outputtrimreads

# run bbduk on all pairs of samples
# sample pairs look like:
# RL41_S1_R1_001.fastq.gz
# RL41_S1_R2_001.fastq.gz
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


    echo "writing output to"
    echo $out1
    echo $out2

    # run bbduk and save ALL the output files    
    $bbduk in1=$in1 in2=$in2 out1=$out1 out2=$out2 minlen=$minlen k=25 mink=8 ktrim=r ref=$adaptors hdist=1 overwrite=f qtrim=rl trimq=$trimq t=$threads bhist=$sampleid"bhist.txt" qhist=$sampleid"qhist.txt" gchist=$sampleid"gchist.txt" aqhist=$sampleid"aqhist.txt" lhist=$sampleid"lhist.txt" > $sampleid"bbduk_log.txt"

done

# run fastqc on all the raw and trimmed data files
find $inputf -name '*.fastq.gz' | xargs fastqc  -o $outputrawqc -t $threads
find $outputtrimreads -name '*.fastq.gz' | xargs fastqc  -o $outputtrimqc -t $threads


# map reads to E. grandis reference
# our trimmed files look like: RL41_S1_R1_001_trimmed.fastq.gz  RL41_S1_R2_001_trimmed.fastq.gz
for in1 in $(find $outputtrimreads -name "*R1_001_trimmed.fastq.gz"); do
    in2=${in1%%R1_001_trimmed.fastq.gz}"R2_001_trimmed.fastq.gz"
    echo ""
    echo "mapping with bbmap"
    echo $in1
    echo $in2

    f1=$(basename $in1)

    bbmapout=$outputbase"bbmap/"
    
    id=${f1%%R1_001_trimmed.fastq.gz}

    path=$bbmapout
    stats=$bbmapout"stats_"$id"trimmed.txt"
    covhist=$bbmapout"covhist_"$id"trimmed.txt"
    ihist=$bbmapout"ihist_"$id"trimmed.txt"
    outsambbmap=$bbmapout$id"trimmed.sam"

    echo "running bbmap"
    echo $bbmapout
    echo $path
    $bbmap in1=$in1 in2=$in2 ref=$ref covstats=$stats covhist=$covhist out=$outsambbmap ihist=$ihist path=$path t=$threads

    echo "converting sams to bams"
    outbambbmap=$bbmapout$id"trimmed"
    samtools view -bS -@ $threads $outsambbmap > $outbambbmap".bam"

    echo "sorting and indexing bams"
    samtools sort -@ $threads  $outbambbmap".bam" $outbambbmap
    samtools index $outbambbmap".bam"

    echo "deleting sams and trimmed reads"
    rm $outsambbmap
    rm $in1
    rm $in2
    
    echo "running qualimap on sorted bams"
    outbamdirbbmap=$bbmapout$id"trimmed_qm"
    $qualimap bamqc -bam $outbambbmap".bam" -outdir $outbamdirbbmap -nt $threads
    
done

echo "running indexcov on sorted bams"
echo "indexcov will run on the following chromosomes"
echo $chroms

# TODO - figure out how to run indexcov on the set of BAMs
#bams=$(echo $bbmapout*.bam)

#for chr in $nuc_chroms; do
#    echo $chr
#    $goleft indexcov -p $chr -c $chr $bams
#done




echo "running multiqc"
multiqc $outputbase -o $outputbase
