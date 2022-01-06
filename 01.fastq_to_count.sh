##########################################################################
# File Name: 01.fastq_to_count.sh
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

echo "multiqc begin"
mkdir ${Result}/02.multiqc
time multiqc ${Result}/01.fastqc/*_fastqc.zip -o ${Result}/02.multiqc
echo "multiqc done"

mkdir ${Result}/03.hisat2
mkdir ${Result}/04.samtools
mkdir ${Result}/05.stringtie
mkdir ${Result}/06.merge
mkdir ${Result}/07.htseq
mkdir ${Result}/08.fpkm

for file in `ls ${Data}/*1.clean.fq.gz`;
do
    rowname=${file##*\/};
    sample=${rowname%%_*};
    echo "${sample} hisat2 begin"
    time hisat2 --dta -t -x ${Ref} -1 ${Data}/${sample}_1.clean.fq.gz -2 ${Data}/${sample}_2.clean.fq.gz -S ${Result}/03.hisat2/${sample}.sam
    echo "${sample} hisat2 done"
    echo "${sample} samtools begin"
    time samtools view -bS ${Result}/03.hisat2/${sample}.sam > ${Result}/04.samtools/${sample}.bam
    time samtools sort ${Result}/04.samtools/${sample}.bam -o ${Result}/04.samtools/Sorted-${sample}.bam
    time samtools sort -n ${Result}/04.samtools/${sample}.bam > ${Result}/04.samtools/Name-Sorted-${sample}.bam
    rm ${Result}/03.hisat2/${sample}.sam
    rm ${Result}/04.samtools/${sample}.bam
    time samtools index ${Result}/04.samtools/Sorted-${sample}.bam
    echo "${sample} samtools done"
    echo "${sample} RSeQC begin"
    time bam_stat.py -i ${Result}/04.samtools/Sorted-${sample}.bam
    time read_distribution.py -i ${Result}/04.samtools/Sorted-${sample}.bam -r $Ref_bed
    echo "${sample} RSeQC done"
    echo "${sample} stringtie begin"
    mkdir ${Result}/05.stringtie/${sample}
    time stringtie ${Result}/04.samtools/Sorted-${sample}.bam -l ${sample} -o ${Result}/05.stringtie/${sample}/${sample}.gtf -p 28 -G $Ref_gff -A ${Result}/05.stringtie/${sample}/${sample}-gene_abund.tab -B -e
    echo "${sample} stringtie done"
    echo "${sample} htseq begin"
    time htseq-count -r pos -f bam ${Result}/04.samtools/Sorted-${sample}.bam $Ref_gtf > ${Result}/07.htseq/pos.${sample}.count
    time htseq-count -r name -f bam ${Result}/04.samtools/Name-Sorted-${sample}.bam $Ref_gtf > ${Result}/07.htseq/name.${sample}.count
    echo "${sample} htseq done"
    echo "${sample} bam to fpkm begin"
    time Rscript /bigdata/wangzhang_guest/chenpeng_project/04_pipeline/05_RNAseq/03.bam_to_fpkm.R --bam ${Result}/04.samtools/Name-Sorted-${sample}.bam --gtf $Ref_gtf --output ${Result}/08.fpkm/${sample}
    echo "${sample} bam to fpkm done"
    echo "${sample} change count header begin"
    new_line="gene_id\t${sample}_counts\t${sample}_fpkm\t${sample}_tpm"
    sed -i "1c${new_line}" ${Result}/08.fpkm/${sample}.count
    echo "${sample} change count header done"
done

mkdir ${Result}/06.merge
echo "merge begin"
for dir in `ls ${Result}/05.stringtie`; do file=${Result}/05.stringtie/${dir}*.gtf; echo "${file}"; done > ${Result}/path.txt
time stringtie --merge -p 10 -G $Ref_gtf -o ${Result}/06.merge/total_merged.gtf -l merge ${Result}/path.txt
echo "merge done"

echo "gffcompare begin"
time gffcompare -r $Ref_gtf -G ${Result}/06.merge/total_merged.gtf -o gffcompare_result
echo "gffcompare done"

echo "Count merge begin"
awk '{print $1}' ${Result}/08.fpkm/${sample}.count > ${Result}/all.count
for file in `ls ${Result}/08.fpkm/*.count`;
do
    cat ${Result}/all.count > ${Result}/log.count
    join ${file} ${Result}/log.count > ${Result}/all.count
done
rm ${Result}/log.count
echo "Count merge done"
