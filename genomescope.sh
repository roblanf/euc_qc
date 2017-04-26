inputf="/disks/dacelo/data/raw_data/tree_EM1/Project_SN7001117R_0083_CKulheim_LBronham_Melaleuca/"
outputf="/disks/dacelo/data/QC/Project_SN7001117R_0083_CKulheim_LBronham_Melaleuca/genomescope/"

threads=20 # how many threads do you want to use?
s='100M' 
kmer_size=21 # 21 suggested in the docs - longer gives issues with seq. error
read_length=100 # length of reads in your input fastq files
kmer_max=5000 # ignore kmers represented more the kmer_max times (avoids issues with e.g. cpDNA)

mkdir $outputf
cd $outputf

# first we join fwd and reverse reads into single .fastq.gz files
for in1 in $(find $inputf -name "*R1_001.fastq.gz"); do
    in2=${in1%%R1_001.fastq.gz}"R2_001.fastq.gz"

    id="$(dirname $in1)"
    id="$(basename $id)"
    readf=$id".fastq.gz"
    cat $in1 $in2 > $readf

done

# here we join biological replicates into files, we need this to get sufficient coverage
# to run genomescope
cat Sample_M1* > M1.fastq.gz
cat Sample_M2* > M2.fastq.gz
cat Sample_M3* > M3.fastq.gz
cat Sample_M4* > M4.fastq.gz
cat Sample_M5* > M5.fastq.gz
cat Sample_M6* > M6.fastq.gz
cat Sample_M7* > M7.fastq.gz
cat Sample_M8* > M8.fastq.gz

rm Sample*

# make one fastq with ALL THE READS (since each sample is only ~30x which is potentially a little low)
cat *.fastq.gz > all.fastq.gz

for readf in $(find *.fastq.gz); do

	# count kmers in jellyfish, then make a histogram. NB: use zcat because jellyfish needs raw fastq, not .gz
	zcat $readf | jellyfish count -C -m $kmer_size -s $s -t $threads -o $outputf$readf'.jf' /dev/fd/0
	jellyfish histo -t $threads $outputf$readf'.jf' > $readf'.histo'

	# run genomescope and put the output in a new directory
	id=`echo "$readf" | cut -d'.' -f1`
	genomescope.R $readf'.histo' $kmer_size $read_length $outputf$id $kmer_max verbose

done

# clean up big files
rm *.fastq.gz
rm *.jf
