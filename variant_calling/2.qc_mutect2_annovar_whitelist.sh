#!/bin/bash -l

#resources:
#$ -l h_rt=12:00:0
#$ -l mem=5G
#$ -l tmpfs=10G
#$ -N mutect2_annovar_whitelist
#$ -t 1-384


module load bwa
module load samtools

module load blic-modules
module load bioconda/2020-11-06

module load python
module load fastqc


#to get into IFC directory
cd /home/regmenu/Scratch/IFC_${1}

#sort and index

#to get into IFC directory
cd /home/regmenu/Scratch/IFC_${1}/S${SGE_TASK_ID}-${1}

pwd
for file in bwa_hg38_*.bam; do
	samtools sort $file -o ${file%.*}.sorted.bam
done
for file2 in bwa_hg38_*.sorted.bam; do
	samtools index ${file2}
done

#depth and qc
for file in bwa_hg38*.sorted.bam; do
        samtools depth -b /home/regmenu/Scratch/mutect_files/amplicon_intervals.bed ${file} > ${file%.*}.samtools.depth.csv
        bedtools coverage -header -d -a /home/regmenu/Scratch/mutect_files/amplicon_intervals.bed -b ${file} > ${file%.*}.bedtools.coverage.csv
done


fastqc -t 8 bwa_hg38*.sorted.bam

multiqc .


#recalibrate base score

gatk --java-options "-Xmx4G" BaseRecalibrator -I bwa_hg38_S${SGE_TASK_ID}-${1}.sorted.bam -R /home/regmenu/Scratch/Hg38/GRCh38.u2af1_fix.fa --known-sites /home/regmenu/Scratch/Homo_sapiens_assembly38.dbsnp138.vcf --known-sites /home/regmenu/Scratch/Homo_sapiens_assembly38.known_indels.vcf.gz --known-sites /home/regmenu/Scratch/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz -O S${SGE_TASK_ID}-${1}_recal_data.table

#apply bqsr

gatk --java-options "-Xmx4G" ApplyBQSR -R /home/regmenu/Scratch/Hg38/GRCh38.u2af1_fix.fa  -I bwa_hg38_*.sorted.bam --bqsr-recal-file S${SGE_TASK_ID}-${1}_recal_data.table -O S${SGE_TASK_ID}-${1}_bqsr.bam


#mutect2 
gatk --java-options "-Xmx4G" Mutect2 -R /home/regmenu/Scratch/Hg38/GRCh38.u2af1_fix.fa -I S${SGE_TASK_ID}-${1}_bqsr.bam --germline-resource /home/regmenu/Scratch/af-only-gnomad.hg38.vcf.gz -L /home/regmenu/Scratch/mutect_files/amplicon_intervals.bed --max-reads-per-alignment-start 0 -O S${SGE_TASK_ID}-${1}_m2.vcf.gz --bam-output S${SGE_TASK_ID}-${1}_m2_bamout.bam --f1r2-tar-gz S${SGE_TASK_ID}-${1}_f1r2.tar.gz


#filter without orientation

gatk --java-options "-Xmx4G" GetPileupSummaries -I S${SGE_TASK_ID}-${1}_bqsr.bam -V /home/regmenu/Scratch/mutect_files/small_exac_common_3.hg38.vcf.gz -L /home/regmenu/Scratch/mutect_files/amplicon_intervals.bed -O S${SGE_TASK_ID}-${1}_getpileup
gatk --java-options "-Xmx4G" CalculateContamination -I S${SGE_TASK_ID}-${1}_getpileup --tumor-segmentation S${SGE_TASK_ID}-${1}_tumor_seg -O S${SGE_TASK_ID}-${1}_contamination
gatk --java-options "-Xmx4G" FilterMutectCalls -R /home/regmenu/Scratch/Hg38/GRCh38.u2af1_fix.fa -V S${SGE_TASK_ID}-${1}_m2.vcf.gz --tumor-segmentation S${SGE_TASK_ID}-${1}_tumor_seg --contamination-table S${SGE_TASK_ID}-${1}_contamination --stats S${SGE_TASK_ID}-${1}_m2.vcf.gz.stats --filtering-stats S${SGE_TASK_ID}-${1}_filt_stats -O S${SGE_TASK_ID}-${1}_m2filtered



gatk --java-options "-Xmx4G" SelectVariants -R /home/regmenu/Scratch/Hg38/GRCh38.u2af1_fix.fa -V S${SGE_TASK_ID}-${1}_m2filtered  -select "DP>50" -O S${SGE_TASK_ID}-${1}_m2filtered2

#annovar

perl /home/regmenu/Scratch/annovar/table_annovar.pl S${SGE_TASK_ID}-${1}_m2filtered2 /home/regmenu/Scratch/annovar/humandb -buildver hg38 -out S${SGE_TASK_ID}-${1}_annovar.vcf -remove -protocol refGene,cosmic70 -operation g,f -nastring . -vcfinput

##whitelist filter

module -f unload compilers mpi gcc-libs
module load r/recommended
export R_LIBS=/home/regmenu/Scratch/MyRlibs_4.2:$R_LIBS

R --vanilla  < /home/regmenu/Scratch/whitelist_filter.R $1 $SGE_TASK_ID

