#!/bin/bash

#SBATCH --job-name=GenoGVCF                   #This is the name of your job
#SBATCH --cpus-per-task=1                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=64G              #This is the memory reserved per core.
#Total memory reserved: 64GB

#SBATCH --time=72:00:00        #This is the time that your task will run
#SBATCH --qos=1week           #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=/scicore/home/schiera/gizevo30/projects/astyanax_var/scripts/logs/GenoGVCFout.txt     #These are the STDOUT and STDERR files
#SBATCH --error=/scicore/home/schiera/gizevo30/projects/astyanax_var/scripts/logs/GenoGVCFerr.txt
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

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T GenotypeGVCFs -R /scicore/home/schiera/gizevo30/projects/astyanax_var/genome/Astyanax_mexicanus.AstMex102.dna.toplevel.fa --variant SRR1575270_variants.g.vcf --variant SRR1575271_variants.g.vcf --variant SRR1575272_variants.g.vcf --variant SRR1575273_variants.g.vcf --variant SRR1575274_variants.g.vcf --variant SRR1575275_variants.g.vcf --variant SRR1575276_variants.g.vcf --variant SRR1575277_variants.g.vcf --variant SRR1575278_variants.g.vcf --variant SRR1575279_variants.g.vcf --variant SRR1575280_variants.g.vcf --variant SRR1575281_variants.g.vcf --variant SRR1575282_variants.g.vcf --variant SRR1575283_variants.g.vcf --variant SRR1575284_variants.g.vcf --variant SRR1575285_variants.g.vcf --variant SRR1575286_variants.g.vcf --variant SRR1575287_variants.g.vcf --variant SRR1575288_variants.g.vcf --variant SRR1575289_variants.g.vcf --variant SRR1575290_variants.g.vcf --variant SRR1575291_variants.g.vcf --variant SRR1575292_variants.g.vcf --variant SRR1575293_variants.g.vcf --variant SRR1575294_variants.g.vcf --variant SRR1575295_variants.g.vcf --variant SRR1575296_variants.g.vcf --variant SRR1575297_variants.g.vcf --variant SRR1575298_variants.g.vcf --variant SRR1575299_variants.g.vcf --variant SRR1927184_variants.g.vcf --variant SRR1927212_variants.g.vcf --variant SRR1927214_variants.g.vcf --variant SRR1927215_variants.g.vcf --variant SRR1927218_variants.g.vcf --variant SRR1927221_variants.g.vcf --variant SRR1927224_variants.g.vcf --variant SRR1927228_variants.g.vcf --variant SRR1927232_variants.g.vcf --variant SRR1927233_variants.g.vcf --variant SRR1927234_variants.g.vcf --variant SRR1927235_variants.g.vcf --variant SRR1927236_variants.g.vcf --variant SRR1927237_variants.g.vcf --variant SRR1927238_variants.g.vcf --max_alternate_alleles 10 -o astyanax_variants.vcf

