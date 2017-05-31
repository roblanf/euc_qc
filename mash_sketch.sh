# this script submits a separate mash sketch job for each reads file 

PROJ="/short/xf1/Epauc"
INDIR=${PROJ}/raw/SN877_merged
OUTDIR=${PROJ}/QC/mash

ls ${INDIR}/*.fastq.gz | while read -r file; do 
	qsub -v INFILE=$file ~/jobs/QC_MASH_SKETCH.job 
done


