
inputf="/disks/dacelo/data/raw_data/tree_EM1/Project_SN7001117R_0083_CKulheim_LBronham_Melaleuca/"
outputf="/disks/dacelo/data/QC/Project_SN7001117R_0083_CKulheim_LBronham_Melaleuca/mash/"

mkdir $outputf

threads=5 # number of threads to use
kmer_size=16 # kmer size to use for minhash sketches
m=2 # Minimum copies of each k-mer required to pass noise filter for reads
s=10000 #Sketch size. Each sketch will have at most this many non-redundant min-hashes

cd $outputf

for in1 in $(find $inputf -name "*R1_001.fastq.gz"); do
    in2=${in1%%R1_001.fastq.gz}"R2_001.fastq.gz"

    id="$(dirname $in1)"
    id="$(basename $id)"
    readf=$id".fastq.gz"
    cat $in1 $in2 > $readf

done

# make the mash sketches
gzs=$(find *.fastq.gz)
mash sketch -p $threads -m $m -k $kmer_size -s $s $gzs

# use mash paste to join the sketches together
fs=$(find *.msh)
mash paste raw $fs   

# verify it worked
mash info raw.msh

# calculate all pairwise distances
mash dist raw.msh raw.msh > distances.tab

Rscript --vanilla heatmap.r $inputf $outputf"heatmap.pdf"