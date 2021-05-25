#!/bin/bash

#SBATCH --job-name=VarFilterSelectVar                   #This is the name of your job
#SBATCH --cpus-per-task=4                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=16G              #This is the memory reserved per core.
#Total memory reserved: 16GB

#SBATCH --time=24:00:00        #This is the time that your task will run
#SBATCH --qos=1day           #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=/scicore/home/schiera/gizevo30/projects/astyanax_var/scripts/logs/VarFilterSelectVarout.txt     #These are the STDOUT and STDERR files
#SBATCH --error=/scicore/home/schiera/gizevo30/projects/astyanax_var/scripts/logs/VarFilterSelectVarerr.txt
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=max.shafer@gmail.com        #You will be notified via email when your task ends or fails

#This job runs from the current working directory


#Remember:
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $JOB_ID stores the ID number of your task.


#load your required modules below
#################################
module load GATK/3.6-Java-1.8.0_92

#export your required environment variables below
#################################################


#add your command lines below
#############################

vcftools --vcf filtered_astyanax_indels.vcf --out vcftools_fst_output.vcf --max-missing 0.5 --weir-fst-pop

