#!/bin/bash

#SBATCH --job-name=StoF_BWAMEM_MergeBam                   #This is the name of your job
#SBATCH --cpus-per-task=4                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=16G              #This is the memory reserved per core.
#Total memory reserved: 64GB

#SBATCH --time=48:00:00        #This is the time that your task will run
#SBATCH --qos=1week           #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=/scicore/home/schiera/gizevo30/projects/astyanax_var/scripts/logs/StoF_BWAMEM_MergeBam_stdout.txt     #These are the STDOUT and STDERR files
#SBATCH --error=/scicore/home/schiera/gizevo30/projects/astyanax_var/scripts/logs/StoF_BWAMEM_MergeBam_stderr.txt

#You selected an array of jobs from 1 to 18 with 18 simultaneous joss
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
module load bwakit/0.7.15_x64-linux

#export your required environment variables below
#################################################


#add your command lines below
#############################

# Step 3

set -o pipefail

java -Xmx8G -jar $EBROOTPICARD/picard.jar SamToFastq I=SRR15752${SLURM_ARRAY_TASK_ID}_marked.bam FASTQ=/dev/stdout CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true TMP_DIR=$TMPDIR | bwa mem -M -t 7 -p /scicore/home/schiera/gizevo30/projects/astyanax_var/genome/Astyanax_mexicanus.AstMex102.dna.toplevel.fa /dev/stdin | java -Xmx16G -jar $EBROOTPICARD/picard.jar MergeBamAlignment ALIGNED_BAM=/dev/stdin UNMAPPED_BAM=SRR15752${SLURM_ARRAY_TASK_ID}_marked.bam OUTPUT=SRR15752${SLURM_ARRAY_TASK_ID}_piped.bam R=/scicore/home/schiera/gizevo30/projects/astyanax_var/genome/Astyanax_mexicanus.AstMex102.dna.toplevel.fa CREATE_INDEX=true ADD_MATE_CIGAR=true CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS TMP_DIR=$TMPDIR

#rm SRR15752${SLURM_ARRAY_TASK_ID}_marked.bam

