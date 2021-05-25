#!/bin/bash

#SBATCH --job-name=SortMarkDups                   #This is the name of your job
#SBATCH --cpus-per-task=1                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=64G              #This is the memory reserved per core.
#Total memory reserved: 64GB

#SBATCH --time=06:00:00        #This is the time that your task will run
#SBATCH --qos=6hours           #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=/scicore/home/schiera/gizevo30/projects/astyanax_var/scripts/logs/SortMarkDups_stdout_70-99.txt     #These are the STDOUT and STDERR files
#SBATCH --error=/scicore/home/schiera/gizevo30/projects/astyanax_var/scripts/logs/SortMarkDups_stderr_70-99.txt

#You selected an array of jobs from 1 to 9 with 9 simultaneous jobs
#SBATCH --array=90%1
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

java -Xmx20G -jar $EBROOTPICARD/picard.jar SortSam INPUT=SRR15752${SLURM_ARRAY_TASK_ID}_piped.bam OUTPUT=SRR15752${SLURM_ARRAY_TASK_ID}_sorted_reads.bam SORT_ORDER=coordinate TMP_DIR=$TMPDIR

java -Xmx20G -jar $EBROOTPICARD/picard.jar MarkDuplicates INPUT=SRR15752${SLURM_ARRAY_TASK_ID}_sorted_reads.bam OUTPUT=SRR15752${SLURM_ARRAY_TASK_ID}_dedup.bam METRICS_FILE=SRR15752${SLURM_ARRAY_TASK_ID}_dedup_metrics.txt TMP_DIR=$TMPDIR

java -Xmx20G -jar $EBROOTPICARD/picard.jar BuildBamIndex INPUT=SRR15752${SLURM_ARRAY_TASK_ID}_dedup.bam
