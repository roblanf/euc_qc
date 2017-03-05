# Basic quality control for Eucalyptus illumina data

# Rob Lanfear, December 2016

# A few things to set before you go

inputf="/disks/dacelo/data/raw_data/Project_SN7001117R_0083_CKulheim_LBronham_Melaleuca/"
outputbase="/disks/dacelo/data/QC/test/"
ref="/disks/dacelo/data/raw_data/active_refs/Emel.fa.gz" # reference file as a fasta
adaptors="/disks/dacelo/data/programs/bbmap/resources/adapters.fa"

# for more control of the parameters, you can change them in the commandlines below
# of course, if there's a parameter you really want to explore, make it a variable here.

threads=50 # number of threads to use
minlen=50 # minimum length of read to keep after trimming
trimq=20 # trim bases with quality < this 

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
mkdir $outputbase"ngm"
mkdir $outputrawqc
mkdir $outputtrimqc
mkdir $outputtrimreads

# run bbduk on all pairs of samples
# sample pairs look like:
# RL41_S1_R1_001.fastq.gz
# RL41_S1_R2_001.fastq.gz
echo "Trimming with bbduk"
date

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
    bbduk.sh in1=$in1 in2=$in2 out1=$out1 out2=$out2 minlen=$minlen k=25 mink=8 ktrim=r ref=$adaptors hdist=1 overwrite=f qtrim=rl trimq=$trimq t=$threads bhist=$sampleid"bhist.txt" qhist=$sampleid"qhist.txt" gchist=$sampleid"gchist.txt" aqhist=$sampleid"aqhist.txt" lhist=$sampleid"lhist.txt" > $sampleid"bbduk_log.txt"

done

echo "Done trimming"
date

# run fastqc on all the raw and trimmed data files
echo "Running fastqc"
date
find $inputf -name '*.fastq.gz' | xargs fastqc  -o $outputrawqc -t $threads
find $outputtrimreads -name '*.fastq.gz' | xargs fastqc  -o $outputtrimqc -t $threads

echo "Done running fastqc"
date

# map reads to E. grandis reference
# our trimmed files look like: RL41_S1_R1_001_trimmed.fastq.gz  RL41_S1_R2_001_trimmed.fastq.gz
echo "Mapping"
date

for in1 in $(find $outputtrimreads -name "*R1_001_trimmed.fastq.gz"); do
    in2=${in1%%R1_001_trimmed.fastq.gz}"R2_001_trimmed.fastq.gz"
    echo ""
    echo "mapping with NextGenMap"
    echo $in1
    echo $in2

    f1=$(basename $in1)

    ngmout=$outputbase"ngm/"
    
    id=${f1%%R1_001_trimmed.fastq.gz}

    stats=$ngmout"stats_"$id"trimmed.txt"
    covhist=$ngmout"covhist_"$id"trimmed.txt"
    ihist=$ngmout"ihist_"$id"trimmed.txt"
    outsamngm=$ngmout$id"trimmed.sam"

    echo "running ngm"
    echo $ngmout
    ngm -t $threads -p -r $ref -1 $in1 -2 $in2 -o $outsamngm

    echo "converting sam to bam"
    outbamngm=$ngmout$id"trimmed"
    samtools view -bS -@ $threads $outsamngm > $outbamngm".bam"

    echo "sorting bam"
    samtools sort -@ $threads  $outbamngm".bam" $outbamngm

    echo "deleting sams and trimmed reads"
    rm $outsamngm
    rm $in1
    rm $in2
    
    echo "running qualimap on sorted bams"
    outbamdirngm=$ngmout$id"trimmed_qm"
    qualimap bamqc -bam $outbamngm".bam" -outdir $outbamdirngm -nt $threads
    
done
echo "Done Mapping"
date

echo "indexing bams"
date
ls ${ngmout}*.bam | parallel "samtools index {}"
echo "done indexing bams"
date

echo "running indexcov"
date
goleft indexcov --directory $indexcov --sex "" ${ngmout}"*.bam"
echo "done running indexcov"
date

echo "running multiqc"
date
multiqc $outputbase -o $outputbase
date
echo "done running multiqc"
date
