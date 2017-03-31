inputf="/disks/dacelo/data/raw_data/tree_EM1/Project_SN7001117R_0083_CKulheim_LBronham_Melaleuca/"
outputf="/disks/dacelo/data/QC/Project_SN7001117R_0083_CKulheim_LBronham_Melaleuca/genomescope/"

threads=20 # how many threads do you want to use?
mem=100000000000 #1000000000 = 1GB; 100000000000 = 100GB
kmer_size=21 # 21 suggested in the docs - longer gives issues with seq. error
read_length=100 # length of reads in your input fastq files
kmer_max=1000 # ignore kmers represented more the kmer_max times (avoids issues with e.g. cpDNA)

fastq_files=$(find $inputf -name "*.fastq.gz")

mkdir $outputf
cd $outputf

# get kmer histogram from jellyfish
jellyfish count -C -m $kmer_size -s $mem -t $threads $fastq_files -o reads.jf
jellyfish histo -t $threads reads.jf > reads.histo

# run genomeScope
Rscript genomescope.R reads.histo $kmer_size $read_length $outputf $kmer_max verbose

