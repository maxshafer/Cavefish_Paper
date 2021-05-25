#!/bin/bash

#SBATCH --job-name=HaploCallerInf                   #This is the name of your job
#SBATCH --cpus-per-task=4                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=6G              #This is the memory reserved per core.
#Total memory reserved: 24GB

#SBATCH --time=70:00:00        #This is the time that your task will run
#SBATCH --qos=infinite           #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=/scicore/home/schiera/gizevo30/projects/astyanax_var/scripts/logs/HaploCallerG_%a_stdout.txt     #These are the STDOUT and STDERR files
#SBATCH --error=/scicore/home/schiera/gizevo30/projects/astyanax_var/scripts/logs/HaploCallerG_%a_stderr.txt

#You selected an array of jobs from 1 to 9 with 9 simultaneous jobs
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

module load GATK/3.7-0-Java-1.8.0_92

#export your required environment variables below
#################################################


#add your command lines below
#############################

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R /scicore/home/schiera/gizevo30/projects/astyanax_var/genome/Astyanax_mexicanus.AstMex102.dna.toplevel.fa -I SRR1927${SLURM_ARRAY_TASK_ID}_realign.bam -nct 4 -ERC GVCF -minPruning 1 -minDanglingBranchLength 1 --max_alternate_alleles 10 --genotyping_mode DISCOVERY -o SRR1927${SLURM_ARRAY_TASK_ID}_variants.g.vcf

