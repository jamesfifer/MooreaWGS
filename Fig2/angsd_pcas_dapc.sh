#creating mafs for ngsrelate
#Maybamscl
ID="/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/Intermediate"
src="/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/SRC"
cd $ID

for i in Maybamscl; do WC=$(wc -l $i | cut -f1 -d " "); WC08=$(echo "scale=2; $WC*0.8"|bc); WC08Final=$(echo $WC08|awk '{print int($1-0.5)}'); echo -e "\
#!/bin/bash \n#$ -V # inherit the submission environment \n#$ -cwd # start job in submission directory
#$ -N MayOctBaypass$i  # job name, anything you want
#$ -P davieslab
#$ -l h_rt=48:00:00 #maximum run time
#$ -M james.e.fifer@gmail.com #your email
#$ -m be
#$ -pe omp 10 \n
ID=\"/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/Intermediate\"
cd \$ID
GENOME_FASTA=/projectnb/wgsgrp/Reference/Ahyacinthus_genome_V1/Ahyacinthus.chrsV1.fasta

FILTERS=\"-uniqueOnly 1 -remove_bads 1 -minMapQ 30 -minQ 35 -dosnpstat 1 -doHWE 1 -sb_pval 1e-3 -maxHetFreq 0.5 -hetbias_pval 1e-3 -skipTriallelic 1 -anc \$GENOME_FASTA -ref \$GENOME_FASTA -snp_pval 1e-5 -minMaf 0.05\"
TODO=\"-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doPost 1 -doGlf 3\"

/projectnb/davieslab/jfifer/Japan_rad/angsd/angsd -b $i -GL 2 -P 10 -minInd $WC08Final \$FILTERS \$TODO -out MayOnlyAllLoci4Relate">$src/MayOctJuvAllLoci.sh; done
#
ID="/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/Intermediate"
src="/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/SRC"
cd $ID

for i in redoMayOctJuvbamscl; do WC=$(wc -l $i | cut -f1 -d " "); WC08=$(echo "scale=2; $WC*0.8"|bc); WC08Final=$(echo $WC08|awk '{print int($1-0.5)}'); echo -e "\
#!/bin/bash \n#$ -V # inherit the submission environment \n#$ -cwd # start job in submission directory
#$ -N MayOctBaypass$i  # job name, anything you want
#$ -P davieslab
#$ -l h_rt=48:00:00 #maximum run time
#$ -M james.e.fifer@gmail.com #your email
#$ -m be
#$ -pe omp 10 \n
ID=\"/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/Intermediate\"
cd \$ID
GENOME_FASTA=/projectnb/wgsgrp/Reference/Ahyacinthus_genome_V1/Ahyacinthus.chrsV1.fasta

FILTERS=\"-uniqueOnly 1 -remove_bads 1 -minMapQ 30 -minQ 35 -dosnpstat 1 -doHWE 1 -sb_pval 1e-3 -maxHetFreq 0.5 -hetbias_pval 1e-3 -skipTriallelic 1 -anc \$GENOME_FASTA -ref \$GENOME_FASTA -snp_pval 1e-5 -minMaf 0.05\"
TODO=\"-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doPost 1 -doGlf 3\"

/projectnb/davieslab/jfifer/Japan_rad/angsd/angsd -b $i -GL 2 -P 10 -minInd $WC08Final \$FILTERS \$TODO -out MayOctJuvOnlyAllLoci4Relate">$src/MayOctJuvAllLoci.sh; done

#run ngsRelate
ID="/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/Intermediate"
src="/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/SRC"
cd $ID


for i in Maybamscl; do WC=$(wc -l $i | cut -f1 -d " "); echo -e "\
#!/bin/bash \n#$ -V # inherit the submission environment \n#$ -cwd # start job in submission directory
#$ -N MayOctBaypass$i  # job name, anything you want
#$ -P davieslab
#$ -l h_rt=48:00:00 #maximum run time
#$ -M james.e.fifer@gmail.com #your email
#$ -m be
#$ -pe omp 10 \n
ID=\"/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/Intermediate\"
cd \$ID
zcat MayOnlyAllLoci4Relate.mafs.gz | cut -f7 |sed 1d >MayOnlyAllLoci4Relate.freq
### run NgsRelate
ngsRelate -p 10 -g MayOnlyAllLoci4Relate.glf.gz -n $WC -f MayOnlyAllLoci4Relate.freq  -O MayOnlyAllLoci4Relate.newres" >$src/runngsrelate.sh; done

ID="/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/Intermediate"
src="/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/SRC"
cd $ID

for i in redoMayOctJuvbamscl; do WC=$(wc -l $i | cut -f1 -d " "); echo -e "\
#!/bin/bash \n#$ -V # inherit the submission environment \n#$ -cwd # start job in submission directory
#$ -N MayOctBaypass$i  # job name, anything you want
#$ -P davieslab
#$ -l h_rt=48:00:00 #maximum run time
#$ -M james.e.fifer@gmail.com #your email
#$ -m be
#$ -pe omp 10 \n
ID=\"/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/Intermediate\"
cd \$ID
zcat MayOctJuvOnlyAllLoci4Relate.mafs.gz | cut -f7 |sed 1d >MayOctJuvOnlyAllLoci4Relate.freq
### run NgsRelate
ngsRelate -p 10 -g MayOctJuvOnlyAllLoci4Relate.glf.gz -n $WC -f MayOctJuvOnlyAllLoci4Relate.freq  -O MayOctJuvOnlyAllLoci4Relate.newres" >$src/runngsrelate.sh; done


##Generate a bams file that contains no sibs and then run angsd to generate the covariance matrix used in PCoAs
##Using the no sibs subsets
#redos removing outlier and extra clone (see angsd_ibs_pca_JF.R)
ID="/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/Intermediate"
src="/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/SRC"
cd $ID
for i in *bamscl*sib*; do WC=$(wc -l $i | cut -f1 -d " "); WC08=$(echo "scale=2; $WC*0.8"|bc); WC08Final=$(echo $WC08|awk '{print int($1-0.5)}'); echo -e "\
#!/bin/bash \n#$ -V # inherit the submission environment \n#$ -cwd # start job in submission directory
#$ -N MayOctBaypass$i  # job name, anything you want
#$ -P davieslab
#$ -l h_rt=48:00:00 #maximum run time
#$ -M james.e.fifer@gmail.com #your email
#$ -m be
#$ -pe omp 10 \n
ID=\"/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/Intermediate\"
cd \$ID
GENOME_FASTA=/projectnb/wgsgrp/Reference/Ahyacinthus_genome_V1/Ahyacinthus.chrsV1.fasta

FILTERS=\"-uniqueOnly 1 -remove_bads 1 -minMapQ 30 -minQ 35 -dosnpstat 1 -doHWE 1 -sb_pval 1e-3 -maxHetFreq 0.5 -hetbias_pval 1e-3 -skipTriallelic 1 -anc \$GENOME_FASTA -ref \$GENOME_FASTA -snp_pval 1e-5 -minMaf 0.05\"
TODO=\"-doMajorMinor 1 -dobcf 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doPost 1 -doGlf 2\"

/projectnb/davieslab/jfifer/Japan_rad/angsd/angsd  -b $i -GL 1 -P 10 -minInd $WC08Final \$FILTERS \$TODO -out $i">$src/${i}.sh; done


#############################Figure 2C waterson's theta ##################
src="/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/SRC"
int="/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/Intermediate"

for i in $int/LTER1.forereef*bamscl.txt; do $src/downsizebams.sh $i $i.downsize 19; done
for i in $int/LTER2.forereef*bamscl.txt; do $src/downsizebams.sh $i $i.downsize 14; done
for i in $int/LTER3.forereef*bamscl.txt; do $src/downsizebams.sh $i $i.downsize 19; done
for i in $int/LTER5.forereef*bamscl.txt; do $src/downsizebams.sh $i $i.downsize 16; done

ID="/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/Intermediate"
src="/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/SRC"
cd $ID

for i in MayOctJuvbamscl; do WC=$(wc -l $i | cut -f1 -d " "); WC08=$(echo "scale=2; $WC*0.8"|bc); WC08Final=$(echo $WC08|awk '{print int($1-0.5)}'); echo -e "\
#!/bin/bash \n#$ -V # inherit the submission environment \n#$ -cwd # start job in submission directory
#$ -N MayOctJuvAllLoci_nominmaf  # job name, anything you want
#$ -P davieslab
#$ -l h_rt=48:00:00 #maximum run time
#$ -M james.e.fifer@gmail.com #your email
#$ -m be
#$ -pe omp 28 \n
ID=\"/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/Intermediate\"
cd \$ID
GENOME_FASTA=/projectnb/wgsgrp/Reference/Ahyacinthus_genome_V1/Ahyacinthus.chrsV1.fasta

FILTERS=\"-uniqueOnly 1 -remove_bads 1 -minMapQ 30 -minQ 35 -dosnpstat 1 -doHWE 1 -sb_pval 1e-3 -maxHetFreq 0.5 -hetbias_pval 1e-3 -skipTriallelic 1 -anc \$GENOME_FASTA -ref \$GENOME_FASTA\"
TODO=\"-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doPost 1 -doGlf 2\"

/projectnb/davieslab/jfifer/Japan_rad/angsd/angsd -b $i -GL 1 -P 28 -minInd $WC08Final \$FILTERS \$TODO -out MayOctJuvAllLoci_nominmaf">$src/MayOctJuvAllLoci_nominmaf.sh; done

zcat MayOctJuvAllLoci_nominmaf.geno.gz | cut -f1,2  >MayOctJuvAllLoci_nominmaf.sites
/projectnb/davieslab/jfifer/Japan_rad/angsd/angsd sites index MayOctJuvAllLoci_nominmaf.sites


ID="/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/Intermediate"
src="/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/SRC"
cd $ID

for i in *bamscl.txt.downsize; do WC=$(wc -l $i | cut -f1 -d " "); WC08=$(echo "scale=2; $WC*0.8"|bc); WC08Final=$(echo $WC08|awk '{print int($1-0.5)}'); echo -e "\
#!/bin/bash \n#$ -V # inherit the submission environment \n#$ -cwd # start job in submission directory
#$ -N theta$i  # job name, anything you want
#$ -P davieslab
#$ -l h_rt=12:00:00 #maximum run time
#$ -M james.e.fifer@gmail.com #your email
#$ -m be
#$ -pe omp 10 \n
ID=\"/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/Intermediate\"
cd \$ID
GENOME_FASTA=/projectnb/wgsgrp/Reference/Ahyacinthus_genome_V1/Ahyacinthus.chrsV1.fasta
TODO=\"-doSaf 1 -anc \$GENOME_FASTA -ref \$GENOME_FASTA\"
/projectnb/davieslab/jfifer/Japan_rad/angsd/angsd -b $i -GL 1 -sites MayOctJuvAllLoci_nominmaf.sites -P 10 \$TODO -out ${i/%bamscl.txt.downsize}.4theta">$src/theta$i.sh; done


ID="/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/Intermediate/"
src="/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/SRC"
cd $ID

>$src/runwinsfs
for i in *4theta.saf.idx; do echo "winsfs -vv $i > ${i/%.saf.idx/}.sfs">>$src/runwinsfs; done

WC=$(wc -l < "$src/runwinsfs")
echo -e "
#!/bin/bash \n#$ -V # inherit the submission environment \n#$ -cwd # start job in submission directory
#$ -N makesfs  # job name, anything you want
#$ -P davieslab
#$ -l h_rt=12:00:00 #maximum run time
#$ -t 1-$WC:1

ID=\"/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/Intermediate\"
cd \$ID
readarray -t CMD_ARRAY <  /projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/SRC/runwinsfs

# bash indexes arrays from 0 but SGE_TASK_ID must start from 1.  Subtract 1 from the
# SGE_TASK_ID to get the correct index into CMD_ARRAY
INDEX=\$((\$SGE_TASK_ID-1))

# Execute the command for this job array
echo Executing command:  \"\${CMD_ARRAY[\$INDEX]}\"
eval \"\${CMD_ARRAY[\$INDEX]}\"
">$src/runwinsfs.sh

#fold the sfs
for i in *4theta.sfs; do winsfs view --fold $i >${i/%4theta.sfs}4thetafold.sfs; done


#reformat sfs files to be used for angsd again (gets rid of header)
for i in *4thetafold.sfs; do sed -i -e "1d" $i; done

#theta calc using sfs as a prior

ID="/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/Intermediate"
src="/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/SRC"

cd $ID
for i in *bamscl.txt.downsize; do echo -e "\
#!/bin/bash \n#$ -V # inherit the submission environment \n#$ -cwd # start job in submission directory
#$ -N theta$i  # job name, anything you want
#$ -P davieslab
#$ -l h_rt=12:00:00 #maximum run time
#$ -M james.e.fifer@gmail.com #your email
#$ -m be
#$ -pe omp 10 \n
ID=\"/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/Intermediate\"
cd \$ID
/projectnb/davieslab/jfifer/Japan_rad/angsd/misc/realSFS saf2theta ${i/%bamscl.txt.downsize}.4theta.saf.idx -sfs ${i/%bamscl.txt.downsize}.4thetafold.sfs -outname ${i/%bamscl.txt.downsize}.4thetafold -P 10 ">$src/theta$i.sh; done

for i in *thetas.idx; do /projectnb/davieslab/jfifer/Japan_rad/angsd/misc/thetaStat do_stat $i; done

