# FastRNAseq
### This is a RNAseq analyse pipeline
> 01.fastq_to_count.sh

This script is a shell script. 

| Software | function | Version |
| ------------ | ------------ | ----------- |
| FastQC | | v0.11.9 |
| MultiQC | | v1.10 |
| hisat2 | | v4.8.2 |
| samtools | | v1.9 |
| RSeQC | bam_stat.py | v4.0.0 |
|  | read_distribution.py | v4.0.0 |
| stringtie | | v2.1.5 |
| htseq-count |  | v0.13.5 |
| gffcompare |  | v0.12.2 |

| R || v4.0.3 |
| SoupX || v1.5.2 |
| DoubletFinder || v2.0.3 |
> 02.count.R
