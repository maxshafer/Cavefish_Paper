#!/bin/bash

#SBATCH --job-name=VarFilterSelectVar                   #This is the name of your job
#SBATCH --cpus-per-task=1                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=24G              #This is the memory reserved per core.
#Total memory reserved: 24GB

#SBATCH --time=06:00:00        #This is the time that your task will run
#SBATCH --qos=6hours           #You will run in this queue

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

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T SelectVariants -R /scicore/home/schiera/gizevo30/projects/astyanax_var/genome/Astyanax_mexicanus.AstMex102.dna.toplevel.fa -V astyanax_variants.vcf -selectType SNP -o raw_astyanax_snps.vcf

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T VariantFiltration -R /scicore/home/schiera/gizevo30/projects/astyanax_var/genome/Astyanax_mexicanus.AstMex102.dna.toplevel.fa -V raw_astyanax_snps.vcf --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filterName "Herman_snp_filter" -o filtered_astyanax_snps.vcf

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T SelectVariants -R /scicore/home/schiera/gizevo30/projects/astyanax_var/genome/Astyanax_mexicanus.AstMex102.dna.toplevel.fa -V astyanax_variants.vcf -selectType INDEL -o raw_astyanax_indels.vcf

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T VariantFiltration -R /scicore/home/schiera/gizevo30/projects/astyanax_var/genome/Astyanax_mexicanus.AstMex102.dna.toplevel.fa -V raw_astyanax_indels.vcf --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" --filterName "Herman_indel_filter" -o filtered_astyanax_indels.vcf

