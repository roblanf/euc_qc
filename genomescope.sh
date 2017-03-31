inputf="/disks/dacelo/data/raw_data/tree_EM1/Project_SN7001117R_0083_CKulheim_LBronham_Melaleuca/Sample_M1a_index2/"
outputf="/disks/dacelo/data/QC/Project_SN7001117R_0083_CKulheim_LBronham_Melaleuca/genomescope/"

threads=20 # how many threads do you want to use?
s='100M' 
kmer_size=21 # 21 suggested in the docs - longer gives issues with seq. error
read_length=100 # length of reads in your input fastq files
kmer_max=1000 # ignore kmers represented more the kmer_max times (avoids issues with e.g. cpDNA)

mkdir $outputf
cd $outputf

for in1 in $(find $inputf -name "*R1_001.fastq.gz"); do
    in2=${in1%%R1_001.fastq.gz}"R2_001.fastq.gz"

    id="$(dirname $in1)"
    id="$(basename $id)"
    readf=$id".fastq.gz"
    cat $in1 $in2 > $readf

	jellyfish count -C -m $kmer_size -s $s -t $threads $readf -o $outputf$id'reads.jf'
	jellyfish histo -t $threads $id'reads.jf' > $id'reads.histo'

	Rscript genomescope.R $id'reads.histo' $kmer_size $read_length $outputf$id $kmer_max verbose

done