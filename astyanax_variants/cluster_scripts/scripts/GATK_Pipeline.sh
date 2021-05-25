#!/bin/bash

#SBATCH --job-name=GATK_pipline                   #This is the name of your job
#SBATCH --cpus-per-task=7                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=16GG              #This is the memory reserved per core.
#Total memory reserved: 16GB

#SBATCH --time=24:00:00        #This is the time that your task will run
#SBATCH --qos=1day           #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=/scicore/home/schiera/gizevo30/projects/astyanax_var/scripts/logs/GATK_pipeline_stdout.txt     #These are the STDOUT and STDERR files
#SBATCH --error=/scicore/home/schiera/gizevo30/projects/astyanax_var/scripts/logs/GATK_pipeline_stderr.txt

#You selected an array of jobs from 1 to 9 with 9 simultaneous jobs
#SBATCH --array=70-87%18
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=max.shafer@gmail.com        #You will be notified via email when your task ends or fails

#This job runs from the current working directory


#Remember:
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $JOB_ID stores the ID number of your task.


#load your required modules below
#################################

#module load picard
#module load bwakit/0.7.15_x64-linux
module load GATK/3.7-0-Java-1.8.0_92

#export your required environment variables below
#################################################


#add your command lines below
#############################

# Step 1 
#java -Xmx8G -jar $EBROOTPICARD/picard.jar FastqToSam FASTQ=SRR15752${SLURM_ARRAY_TASK_ID}_1.fastq FASTQ2=SRR15752${SLURM_ARRAY_TASK_ID}_2.fastq OUTPUT=SRR15752${SLURM_ARRAY_TASK_ID}_temp.bam SAMPLE_NAME=SRR15752${SLURM_ARRAY_TASK_ID} LIBRARY_NAME=SRR15752${SLURM_ARRAY_TASK_ID} PLATFORM=illumina

# Step 2

#java -Xmx8G -jar $EBROOTPICARD/picard.jar MarkIlluminaAdapters I=SRR15752${SLURM_ARRAY_TASK_ID}_temp.bam O=SRR15752${SLURM_ARRAY_TASK_ID}_marked.bam M=SRR15752${SLURM_ARRAY_TASK_ID}_markilluminaadapters_metrics.txt

#rm SRR15752${SLURM_ARRAY_TASK_ID}_temp.bam
#rm SRR15752${SLURM_ARRAY_TASK_ID}_1.fastq
#rm SRR15752${SLURM_ARRAY_TASK_ID}_2.fastq

# Step 3

#set -o pipefail

#java -Xmx8G -jar $EBROOTPICARD/picard.jar SamToFastq I=SRR15752${SLURM_ARRAY_TASK_ID}_marked.bam FASTQ=/dev/stdout CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true TMP_DIR=/scicore/home/schiera/gizevo30/mstemp | bwa mem -M -t 7 -p /scicore/home/schiera/gizevo30/projects/astyanax_var/genome/Astyanax_mexicanus.AstMex102.dna.toplevel.fa /dev/stdin | java -Xmx16G -jar $EBROOTPICARD/picard.jar MergeBamAlignment ALIGNED_BAM=/dev/stdin UNMAPPED_BAM=SRR15752${SLURM_ARRAY_TASK_ID}_marked.bam OUTPUT=SRR15752${SLURM_ARRAY_TASK_ID}_piped.bam R=/scicore/home/schiera/gizevo30/projects/astyanax_var/genome/Astyanax_mexicanus.AstMex102.dna.toplevel.fa CREATE_INDEX=true ADD_MATE_CIGAR=true CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS TMP_DIR=/scicore/home/schiera/gizevo30/mstemp

#rm SRR15752${SLURM_ARRAY_TASK_ID}_marked.bam

# Step 4

#java -Xmx16G -jar $EBROOTPICARD/picard.jar SortSam INPUT=SRR15752${SLURM_ARRAY_TASK_ID}_piped.bam OUTPUT=SRR15752${SLURM_ARRAY_TASK_ID}_sorted_reads.bam SORT_ORDER=coordinate

#java -Xmx16G -jar $EBROOTPICARD/picard.jar MarkDuplicates INPUT=SRR15752${SLURM_ARRAY_TASK_ID}_sorted_reads.bam OUTPUT=SRR15752${SLURM_ARRAY_TASK_ID}_dedup.bam METRICS_FILE=SRR15752${SLURM_ARRAY_TASK_ID}_dedup_metrics.txt

#java -Xmx16G -jar $EBROOTPICARD/picard.jar BuildBamIndex INPUT=SRR15752${SLURM_ARRAY_TASK_ID}_dedup.bam

rm SRR15752${SLURM_ARRAY_TASK_ID}_piped.bam
rm SRR15752${SLURM_ARRAY_TASK_ID}_sorted_reads.bam

# Step 5

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R /scicore/home/schiera/gizevo30/projects/astyanax_var/genome/Astyanax_mexicanus.AstMex102.dna.toplevel.fa -I SRR15752${SLURM_ARRAY_TASK_ID}_dedup.bam -o SRR15752${SLURM_ARRAY_TASK_ID}_realign.intervals

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T IndelRealigner -R /scicore/home/schiera/gizevo30/projects/astyanax_var/genome/Astyanax_mexicanus.AstMex102.dna.toplevel.fa -targetIntervals SRR15752${SLURM_ARRAY_TASK_ID}_realign.intervals -I SRR15752${SLURM_ARRAY_TASK_ID}_dedup.bam -o SRR15752${SLURM_ARRAY_TASK_ID}_realign.bam

#rm SRR15752${SLURM_ARRAY_TASK_ID}_dedup.bam

# Step 7

#java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R /scicore/home/schiera/gizevo30/projects/astyanax_var/genome/Astyanax_mexicanus.AstMex102.dna.toplevel.fa -I SRR15752${SLURM_ARRAY_TASK_ID}_realign.bam -minPruning 1 -minDanglingBranchLength 1 --max_alternate_alleles 10 --genotyping_mode DISCOVERY -o SRR15752${SLURM_ARRAY_TASK_ID}_variants.g.vcf

