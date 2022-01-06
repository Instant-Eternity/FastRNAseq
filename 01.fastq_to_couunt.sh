##########################################################################
# File Name: 01.fastq_to_couunt.sh
# Author: Instant-Eternity
# mail: hunterfirstone@i.smu.edu.cn
# Created Time: Thu 03 Jun 2021 11:34:14 AM CST
#########################################################################
#!/use/bin/bash
Work_Path=/bigdata/wangzhang_guest/chenpeng_project
Data=${Work_Path}/01_data/19_LJX_RNAseq
mkdir ${Work_Path}/06_result/19_LJX_RNAseq
Result=${Work_Path}/06_result/19_LJX_RNAseq
Ref=${Work_Path}/02_reference/02_ref_mm39/mm39/mm39
Ref_gff=${Work_Path}/02_reference/02_ref_mm39/gencode.vM28.annotation.gff3
Ref_gtf=${Work_Path}/02_reference/02_ref_mm39/gencode.vM28.annotation.gtf
Ref_bed=${Work_Path}/02_reference/02_ref_mm39/mm39.bed

starttime=$(date +%Y-%m-%d\ %H:%M:%S)
echo ${starttime}

time conda activate scRNAseq

echo "fastqc begin"
mkdir ${Result}/01.fastqc
time fastqc -o ${Result}/01.fastqc -f fastq ${Data}/*.fq.gz
echo "fastqc done"
