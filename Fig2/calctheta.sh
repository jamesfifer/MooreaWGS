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


##########

ID="/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/Intermediate"
src="/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/SRC"
cd $ID

cat *May*bamscl.txt.downsize >May.FR.bamscl.downsize
cat *Oct*bamscl.txt.downsize >Oct.FR.bamscl.downsize
cat *NA*bamscl.txt.downsize >NA.FR.bamscl.downsize


ID="/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/Intermediate"
src="/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/SRC"
cd $ID


for i in *FR.bamscl.downsize; do WC=$(wc -l $i | cut -f1 -d " "); WC08=$(echo "scale=2; $WC*0.8"|bc); WC08Final=$(echo $WC08|awk '{print int($1-0.5)}'); echo -e "\
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
/projectnb/davieslab/jfifer/Japan_rad/angsd/angsd -b $i -GL 1 -sites MayOctJuvAllLoci_nominmaf.sites -P 10 \$TODO -out ${i/%FR.bamscl.downsize}.4theta">$src/theta$i.sh; done


for i in theta*FR*sh; do qsub $i; done

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


#reformat sfs files to be used for angsd again (gets rid of header)
for i in *..4theta.sfs; do sed -i -e "1d" $i; done

#theta calc using sfs as a prior

ID="/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/Intermediate"
src="/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/SRC"

cd $ID
for i in *FR.bamscl.downsize; do echo -e "\
#!/bin/bash \n#$ -V # inherit the submission environment \n#$ -cwd # start job in submission directory
#$ -N theta$i  # job name, anything you want
#$ -P davieslab
#$ -l h_rt=12:00:00 #maximum run time
#$ -M james.e.fifer@gmail.com #your email
#$ -m be
#$ -pe omp 10 \n
ID=\"/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/Intermediate\"
cd \$ID
/projectnb/davieslab/jfifer/Japan_rad/angsd/misc/realSFS saf2theta ${i/%FR.bamscl.downsize}.4theta.saf.idx -sfs ${i/%FR.bamscl.downsize}.4theta.sfs -outname ${i/%FR.bamscl.downsize}.4theta -P 10 ">$src/theta$i.sh; done

for i in *thetas.idx; do /projectnb/davieslab/jfifer/Japan_rad/angsd/misc/thetaStat do_stat $i; done




###############
ID="/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/Intermediate"
src="/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/SRC"
cd $ID



for i in MayOctbamscl; do WC=$(wc -l $i | cut -f1 -d " "); WC08=$(echo "scale=2; $WC*0.8"|bc); WC08Final=$(echo $WC08|awk '{print int($1-0.5)}'); echo -e "\
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

/projectnb/davieslab/jfifer/Japan_rad/angsd/angsd -b $i -GL 1 -P 28 -minInd $WC08Final \$FILTERS \$TODO -out MayOctAllLoci_nominmaf">$src/MayOctAllLoci_nominmaf.sh; done

zcat MayOctAllLoci_nominmaf.geno.gz | cut -f1,2  >MayOctAllLoci_nominmaf.sites
/projectnb/davieslab/jfifer/Japan_rad/angsd/angsd sites index MayOctAllLoci_nominmaf.sites


src="/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/SRC"
int="/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/Intermediate"

for i in $int/LTER2.crest*bamscl.txt; do $src/downsizebams.sh $i $i.downsize 13; done

ID="/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/Intermediate"
src="/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/SRC"
cd $ID


while read i; do WC=$(wc -l $i | cut -f1 -d " "); WC08=$(echo "scale=2; $WC*0.8"|bc); WC08Final=$(echo $WC08|awk '{print int($1-0.5)}'); echo -e "\
#!/bin/bash \n#$ -V # inherit the submission environment \n#$ -cwd # start job in submission directory
#$ -N theta$i  # job name, anything you wantLTER1.backreef.Maybamscl.txt.4theta.saf.idxLTER1.backreef.Maybamscl.txt.4theta.saf.idx
#$ -P davieslab
#$ -l h_rt=12:00:00 #maximum run time
#$ -M james.e.fifer@gmail.com #your email
#$ -m be
#$ -pe omp 10 \n
ID=\"/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/Intermediate\"
cd \$ID
GENOME_FASTA=/projectnb/wgsgrp/Reference/Ahyacinthus_genome_V1/Ahyacinthus.chrsV1.fasta
TODO=\"-doSaf 1 -anc \$GENOME_FASTA -ref \$GENOME_FASTA\"
/projectnb/davieslab/jfifer/Japan_rad/angsd/angsd -b $i -GL 1 -sites MayOctJuvAllLoci_nominmaf.sites -P 10 \$TODO -out ${i/%bamscl.txt.downsize}.4theta">$src/theta$i.sh; done <$ID/bams2calctheta


ID="/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/Intermediate/"
src="/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/SRC"
cd $ID

>$src/runwinsfs
while read i; do echo "winsfs -vv ${i/%bamscl.txt.downsize}.4theta.saf.idx > ${i/%bamscl.txt.downsize}.4theta.sfs">>$src/runwinsfs; done <$ID/bams2calctheta


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

while read i; do sed -i -e "1d" ${i/%bamscl.txt.downsize}.4theta.sfs; done <bams2calctheta

ID="/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/Intermediate"
src="/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/SRC"

cd $ID
while read i in; do echo -e "\
#!/bin/bash \n#$ -V # inherit the submission environment \n#$ -cwd # start job in submission directory
#$ -N theta$i  # job name, anything you want
#$ -P davieslab
#$ -l h_rt=12:00:00 #maximum run time
#$ -M james.e.fifer@gmail.com #your email
#$ -m be
#$ -pe omp 10 \n
ID=\"/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/Intermediate\"
cd \$ID
/projectnb/davieslab/jfifer/Japan_rad/angsd/misc/realSFS saf2theta ${i/%bamscl.txt.downsize}.4theta.saf.idx -sfs ${i/%bamscl.txt.downsize}.4theta.sfs -outname ${i/%bamscl.txt.downsize}.4theta -P 10 ">$src/theta$i.sh; done <bams2calctheta

for i in *thetas.idx; do /projectnb/davieslab/jfifer/Japan_rad/angsd/misc/thetaStat do_stat $i; done


