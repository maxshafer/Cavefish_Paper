#!/bin/bash

#SBATCH --job-name=uBAM_MarkAdapt                   #This is the name of your job
#SBATCH --cpus-per-task=1                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=24G              #This is the memory reserved per core.
#Total memory reserved: 24GB

#SBATCH --time=06:00:00        #This is the time that your task will run
#SBATCH --qos=6hours           #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=/scicore/home/schiera/gizevo30/projects/astyanax_var/scripts/logs/uBAM_MarkAdapt_stdout.txt     #These are the STDOUT and STDERR files
#SBATCH --error=/scicore/home/schiera/gizevo30/projects/astyanax_var/scripts/logs/uBAM_MarkAdapt_stderr.txt

#You selected an array of jobs from 1 to 9 with 9 simultaneous jobs
#SBATCH --array=70-99%30
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

java -Xmx8G -jar $EBROOTPICARD/picard.jar FastqToSam FASTQ=SRR15752${SLURM_ARRAY_TASK_ID}_1.fastq FASTQ2=SRR15752${SLURM_ARRAY_TASK_ID}_2.fastq OUTPUT=SRR15752${SLURM_ARRAY_TASK_ID}_temp.bam SAMPLE_NAME=SRR15752${SLURM_ARRAY_TASK_ID} LIBRARY_NAME=SRR15752${SLURM_ARRAY_TASK_ID} PLATFORM=illumina

rm SRR15752${SLURM_ARRAY_TASK_ID}_1.fastq
rm SRR15752${SLURM_ARRAY_TASK_ID}_2.fastq

java -Xmx8G -jar $EBROOTPICARD/picard.jar MarkIlluminaAdapters I=SRR15752${SLURM_ARRAY_TASK_ID}_temp.bam O=SRR15752${SLURM_ARRAY_TASK_ID}_marked.bam M=SRR15752${SLURM_ARRAY_TASK_ID}_markilluminaadapters_metrics.txt

rm SRR15752${SLURM_ARRAY_TASK_ID}_temp.bam

