#!/bin/bash

#SBATCH --job-name=uBAM_MarkAdapt                   #This is the name of your job
#SBATCH --cpus-per-task=1                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=24G              #This is the memory reserved per core.
#Total memory reserved: 24GB

#SBATCH --time=06:00:00        #This is the time that your task will run
#SBATCH --qos=1day           #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=/scicore/home/schiera/gizevo30/projects/astyanax_var/scripts/logs/uBAM_MarkAdapt_stdout.txt     #These are the STDOUT and STDERR files
#SBATCH --error=/scicore/home/schiera/gizevo30/projects/astyanax_var/scripts/logs/uBAM_MarkAdapt_stderr.txt

#You selected an array of jobs from 1 to 18 with 18 simultaneous jobs
#SBATCH --array=184-238%55
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=max.shafer@gmail.com        #You will be notified via email when your task ends or fails

#This job runs from the current working directory


#Remember:
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $JOB_ID stores the ID number of your task.


#load your required modules below
#################################

module load picard

#export your required environment variables below
#################################################


#add your command lines below
#############################

# Step 1 
java -Xmx8G -jar $EBROOTPICARD/picard.jar FastqToSam FASTQ=SRR1927${SLURM_ARRAY_TASK_ID}_1.fastq FASTQ2=SRR1927${SLURM_ARRAY_TASK_ID}_2.fastq OUTPUT=SRR1927${SLURM_ARRAY_TASK_ID}_temp.bam SAMPLE_NAME=SRR1927${SLURM_ARRAY_TASK_ID} LIBRARY_NAME=SRR1927${SLURM_ARRAY_TASK_ID} PLATFORM=illumina TMP_DIR=$TMPDIR

rm SRR1927${SLURM_ARRAY_TASK_ID}_1.fastq
rm SRR1927${SLURM_ARRAY_TASK_ID}_2.fastq

# Step 2

java -Xmx8G -jar $EBROOTPICARD/picard.jar MarkIlluminaAdapters I=SRR1927${SLURM_ARRAY_TASK_ID}_temp.bam O=SRR1927${SLURM_ARRAY_TASK_ID}_marked.bam M=SRR1927${SLURM_ARRAY_TASK_ID}_markilluminaadapters_metrics.txt TMP_DIR=$TMPDIR

rm SRR1927${SLURM_ARRAY_TASK_ID}_temp.bam

