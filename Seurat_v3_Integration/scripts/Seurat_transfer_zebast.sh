#!/bin/bash

#SBATCH --job-name=SeuratInt                   #This is the name of your job
#SBATCH --cpus-per-task=1                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=64G              #This is the memory reserved per core.
#Total memory reserved: 64GB

#SBATCH --time=24:00:00        #This is the time that your task will run
#SBATCH --qos=1day           #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=/scicore/home/schiera/gizevo30/projects/Seurat_v3_Integration/logs/SeuratIntout.txt     #These are the STDOUT and STDERR files
#SBATCH --error=/scicore/home/schiera/gizevo30/projects/Seurat_v3_Integration/logs/SeuratInterr.txt
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=max.shafer@gmail.com        #You will be notified via email when your task ends or fails

#This job runs from the current working directory


#Remember:
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $JOB_ID stores the ID number of your task.


#load your required modules below
#################################

module load R/3.5.2-foss-2018b

#add your command lines below 
#############################  

Rscript Seurat_transfer_zebast.R

