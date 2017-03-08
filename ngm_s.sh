#!/bin/bash
# Compare various settings for sensitivity in NextGenMap

# Rob Lanfear, March 2016

# forward and reverse input reads
R1="/disks/dacelo/data/raw_data/Project_SN7001117R_0083_CKulheim_LBronham_Melaleuca/Sample_M1a_index2/M1a_index2_CGATGTAT_L003_R1_001.fastq.gz"
R2="/disks/dacelo/data/raw_data/Project_SN7001117R_0083_CKulheim_LBronham_Melaleuca/Sample_M1a_index2/M1a_index2_CGATGTAT_L003_R2_001.fastq.gz"

# directory where you want to do stuff
outputbase="/disks/dacelo/data/QC/Project_SN7001117R_0083_CKulheim_LBronham_Melaleuca/compare_ngm_s/"

# reference genome
ref="/disks/dacelo/data/raw_data/active_refs/Emel.fa.gz"
gff="/disks/dacelo/data/raw_data/active_refs/Egrandis_genes_chr1_to_chr11.gff3"
threads=30

declare -a s=("0.1" "0.2" "0.3" "0.4" "0.5" "0.6" "0.7" "0.8" "0.9")

## now loop through the above array to try out all of the s values you want
for i in "${s[@]}"
do
	echo "running ngm with s set to $i"
	outbam="out_"$i.bam
	
	ngm -s $i -t $threads -p -r $ref -1 $R1 -2 $R2 -o "out_"$i.sam 
	samtools view -bS -@ $threads "out_"$i.sam > $outbam
	samtools sort -@ $threads $outbam -o $outbam
	samtools index $outbam
	qualimap bamqc -bam $outbam -outdir $outngm"qualimap_all_"$i"/" -nt $threads -c
	qualimap bamqc -bam $outbam -outdir $outngm"qualimap_gff_"$i"/" -gff $gff -nt $threads -c

done
