#!/bin/bash -l

#resources:
#$ -l h_rt=48:00:0
#$ -l mem=20G
#$ -l tmpfs=300G
#$ -N align_to_bam

module load bwa
module load samtools


#to get into IFC directory
cd /home/regmenu/Scratch/IFC_${1}

for d in */; do
        cd ./${d}
        gunzip *fq.gz
        line=$(head -n 1 *1.fq)


        #remove first character (@):
        line2="${line:1}"

        #remove characters after space in first line of fastq
        line3=${line2%' '*}

        #split by ":"

        IFS=':' read -ra my_array <<< "$line3"
        ##get sample number
        file=(`ls *1.fq`)
        IFS='_' read -r -a sample <<< "$file"
        a=${d::-1}
	bwa mem -t 8 -R "@RG\tID:${my_array[0]}.${my_array[3]}\tPL:ILLUMINA\tPU:${my_array[2]}.${my_array[3]}\tLB:libraries_2-9\tSM:${a}" /home/regmenu/Scratch/Hg38/GRCh38.u2af1_fix.fa /home/regmenu/Scratch/IFC_${1}/${a}/*1.fq /home/regmenu/Scratch/IFC_${1}/${a}/*2.fq | samtools view -o bwa_hg38_${a}.bam
        cd ..
done



