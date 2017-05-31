input_L1="/short/xf1/Epauc/raw/SN877_0428_RLanfear_RSB_Eucalyptus_gDNA"        # lane 1
input_L2="/short/xf1/Epauc/raw/SN877_0431_RLanfear_RSB_Eucalyptus_gDNA"        # Lane 2
outputf="/short/xf1/Epauc/raw/SN877_merged"	# destination directory

mkdir $outputf

for in1 in $(find $input_L1 -name "*R1_001.fastq.gz"); do

        f1=$(basename $in1)
        f2=${f1%%R1_001.fastq.gz}"R2_001.fastq.gz"
        id=$(echo $f1 | cut -f1,2 -d'_')

        cat ${input_L1}/$f1 ${input_L1}/$f2 ${input_L2}/$f1 ${input_L2}/$f2 > ${outputf}/$id".fastq.gz"
done

