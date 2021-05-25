#!/bin/bash

#SBATCH --job-name=fastq_dump                   #This is the name of your job
#SBATCH --cpus-per-task=1                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=1G              #This is the memory reserved per core.
#Total memory reserved: 1GB

#SBATCH --time=01:00:00        #This is the time that your task will run
#SBATCH --qos=6hours           #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=/scicore/home/schiera/gizevo30/projects/astyanax_var/scripts/logs/fastq_dump_stdout.txt     #These are the STDOUT and STDERR files
#SBATCH --error=/scicore/home/schiera/gizevo30/projects/astyanax_var/scripts/logs/fastq_dump_stderr.txt

#You selected an array of jobs from 1 to 18 with 18 simultaneous jobs
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
module load SRA-Toolkit/2.8.1-3-centos_linux64

#export your required environment variables below
#################################################


#add your command lines below
#############################

fastq-dump --outdir /scicore/home/schiera/gizevo30/projects/astyanax_var/sra_reads_nobackup --split-files /scicore/home/schiera/gizevo30/ncbi/public/sra/SRR15752${SLURM_ARRAY_TASK_ID}.sra


