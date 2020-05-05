#### Overview
This sequencing run has 11 human samples and 13 mouse samples that will have to be split up for RSEM alignments to different references.

First pursue RSEM alignments as before and if poor, troubleshoot contamination with HISAT2.

As before, followed sequence similar to [Sham's BCC10 alignment May 2019](Sham's BCC10 alignment May 2019). Only exceptions really noted below.

```bash
cd /scratch/ss4cf
mkdir reseq
mkdir concat
mkdir fastqc
mkdir raw
mkdir rsem
mkdir trim

#In a new tab:
scp -r /Volumes/JanesLab/RNA-seq\ raw\ data/07012019-Janes_Pool12-137912779/* ss4cf@rivanna.hpc.virginia.edu:/scratch/ss4cf/reseq/raw


#Concatenate

pwd
#/scratch/ss4cf/reseq (path name)

# Get the sample names
SAMPLIST=$(find ./raw/*/*/* | cut -f 5 -d "/" | cut -f 1 -d "_" | sort -u)

for SAMP in $SAMPLIST
do
R1=$(find ./raw/*/*/$SAMP*R1*.fastq.gz | sort)
R2=$(find ./raw/*/*/$SAMP*R2*.fastq.gz | sort)
cat $R1 > ./concat/${SAMP}_R1.fastq.gz
cat $R2 > ./concat/${SAMP}_R2.fastq.gz
done

#I usually go check that these files are actually in my concat folder

cd concat
ls

#Trim adapter sequences
#This is a dry run and produces a list of commands, test one to make sure there are no errors
find ./concat/*fastq.gz | cut -f 1 -d "_" | cut -f 3 -d "/" | sort | parallel -j 6 --dry-run "~/software/ea-utils/ea-utils-master/clipper/fastq-mcf -q 10 -t 0.01 -k 0 ~/adapters/nextera_adapters.fa ./concat/{}_R1.fastq.gz ./concat/{}_R2.fastq.gz -o ./trim/{}_R1.trim.fq.gz -o ./trim/{}_R2.trim.fq.gz > ./trim/{}.trim.log"

# Trim adapter sequences - previously did this on the front end, but going to try to do it as a job

#make a trim.sh file:
nano trim.sh

--------
#!/bin/sh

find ./concat/*fastq.gz | cut -f 1 -d "_" | cut -f 3 -d "/" | sort | parallel -j 6 "~/software/ea-utils/ea-utils-master/clipper/fastq-mcf -q 10 -t 0.01 -k 0 ~/adapters/nextera_adapters.fa ./concat/{}_R1.fastq.gz ./concat/{}_R2.fastq.gz -o ./trim/{}_R1.trim.fq.gz -o ./trim/{}_R2.trim.fq.gz > ./trim/{}.trim.log"

---------
#zipped files are necessary for bowtie

sbatch -A janeslab -p standard -t 1:00:00 -n 1 --cpus-per-task 5 ./trim.sh
#can request up to 20, but then your job could wait pending resources etc...


#QC: Check number of reads before and after trimming:

cd trim
SAMPLIST=$(find *log | grep -v _R | cut -f 1 -d "." | sort -u)
for SAMP in $SAMPLIST
do
RC=$(grep "Total reads" ${SAMP}.trim.log | cut -f 3 -d " ")
TR=$(grep "Too short" ${SAMP}.trim.log | cut -f 5 -d " ")
RL=$(( $RC - $TR ))
echo -e '|'$SAMP'|'$RC'|'$RL'|'
done


#Run FASTQC (might have to submit this as separate jobs)
cd ../trim
module load fastqc
#Run FastQC in parallel
find *fq.gz | parallel -j 6 "fastqc {}"
mv ./*fastqc* ../fastqc

```

#### Alignments with RSEM

Will have to separate out human and mouse ones for this, so maybe make two separate align folders and go from there

Fewer CPU requested because it makes jobs get queued for longer than they would have taken anyway!

```bash
cd ..
pwd
#/scratch/ss4cf/reseq (path name)

rm -r rsem
mkdir rsem_h
mkdir rsem_m

cd rsem_h
nano submitjobs.sh
-------------
#!/bin/bash

SAMPLIST=$(find ../trim/BCC*.fq.gz |cut -f 3 -d "/"|  cut -f 1 -d "_" | sort -u)

for SAMP in $SAMPLIST
do

sbatch -A janeslab -p standard -t 10:00:00 -n 1 --cpus-per-task 5 --wrap="sh ./calculate_rsem.sh $SAMP"
sleep 1 # wait 1 second between each job submission

done
------------

nano calculate_rsem.sh
------------
#!/bin/bash/

SAMP=$1

REF="/home/ss4cf/genomes/grch38_cdna_84_ercc/grch38_84_ercc_rsem"

module load rsem

mkdir $SAMP

echo "Processing file $SAMP, starting at $(date +%F\ %T)"

rsem-calculate-expression -p 20 --bowtie2 --bowtie2-path /home/ss4cf/software/bowtie2/bowtie2-2.3.4.3 --single-cell-prior --paired-end ../trim/${SAMP}_R1.trim.fq.gz ../trim/${SAMP}_R2.trim.fq.gz $REF ./${SAMP}/$SAMP

echo "Finished processing file $SAMP at $(date +%F\ %T)"
---------------

chmod u+x *.sh #convert to executable
./submitjobs.sh #run
squeue -u ss4cf #check job queue/status

cd ../rsem_m
nano submitjobs_m.sh
-------------
#!/bin/bash

SAMPLIST=$(find ../trim/*.fq.gz |cut -f 3 -d "/"|  cut -f 1 -d "_" | sort -u | grep -v 'BCC*')

for SAMP in $SAMPLIST
do

sbatch -A janeslab -p standard -t 10:00:00 -n 1 --cpus-per-task 5 --wrap="sh ./calculate_rsem_m.sh $SAMP"
sleep 1 # wait 1 second between each job submission

done
------------

nano calculate_rsem_m.sh
------------
#!/bin/bash/

SAMP=$1

REF="/home/ss4cf/genomes/grcm38_cdna_82_ercc/grcm38ercc_rsem"

module load rsem

mkdir $SAMP

echo "Processing file $SAMP, starting at $(date +%F\ %T)"

rsem-calculate-expression -p 20 --bowtie2 --bowtie2-path /home/ss4cf/software/bowtie2/bowtie2-2.3.4.3 --single-cell-prior --paired-end ../trim/${SAMP}_R1.trim.fq.gz ../trim/${SAMP}_R2.trim.fq.gz $REF ./${SAMP}/$SAMP

echo "Finished processing file $SAMP at $(date +%F\ %T)"
---------------

chmod u+x *.sh #convert to executable
./submitjobs_m.sh #run
squeue -u ss4cf #check job queue/status

#Consolidate alignments before running QC?
cp -r rsem_h/* rsem_m

cd rsem_m
ls

#all RSEM output for 24 samples should be here now


# Load anaconda to use multiqc package which needs some python
module load anaconda
module load bioconda/py3.6
pip install --user multiqc

#Copy trim QC files to RSEM folder to consolidate
cp -r ../fastqc .
multiqc -f . #run multiQC

logout 

#Navigate to folder on your computer where you want to save it
cd /Users/ss4cf/Desktop/Lab/Dylan\ Schaff/RNAseq/KP1_impos 

#Copy multiQC results
scp -r ss4cf@rivanna.hpc.virginia.edu:/scratch/ss4cf/reseq/rsem_m/multiqc* . 
#This includes a nice .html page to view results

```
Nothing was horrendously low so I'm going forward with getting the read counts done for analysis!

```bash
logout
cd /Users/ss4cf/Desktop/Lab/Data/BCC\ Samples/BCC_RNAseq/RSEM\ data/BCC10_KP1ImPos_OPC_reseq/BCC10_KP1ImPos_OPC_reseq_alignments_results 

scp ss4cf@rivanna.hpc.virginia.edu:/scratch/ss4cf/reseq/rsem_m/*/*.isoforms.results .

#I had moved everything to the rsem_m folder to do the QC

#Apparently, KP1ImPos-Pool25 didn't have enough resources to create the .isoforms.results file - so run again.

#submitjob:
SAMPLIST=$(find ../trim/KP1ImPos-Pool25*.fq.gz |cut -f 3 -d "/"|  cut -f 1 -d "_" | sort -u)

sbatch -A janeslab -p standard -t 10:00:00 -n 1 --cpus-per-task 5 --wrap="sh ./calculate_rsem_m.sh $SAMP"

#also change first line of calculate_rsem_m.sh

```
Next use R script for producing RSEM output csv!
