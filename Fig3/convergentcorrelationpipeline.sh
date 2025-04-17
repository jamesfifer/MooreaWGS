#!/bin/bash
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N runpgs  # job name, anything you want
#$ -l h_rt=36:00:00
#$ -M james.e.fifer@gmail.com #your email
#$ -m be
#$ -pe omp 18

module load htslib/1.16
module load bcftools/1.9

module load plink2
bcftools reheader -s new_samples.txt -o MayOctAllLoci_reheader.vcf /projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/Intermediate/MayOctAllLoci_annotated.vcf

plink2 --allow-extra-chr --vcf /projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/Intermediate/MayOctAllLoci_reheader.vcf --double-id --freq --pheno allpops.txt --loop-cats population
