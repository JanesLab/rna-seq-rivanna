Based on `KP1 liver mets` alignments

#FIRST: only ran 1 cartridge to make sure we're getting good alignment rates for all samples, and then the same pool will be run over another cartridge

Use reference genome already prepared for RSEM

#### Set up raw data

```bash
cd /scratch/ss4cf
mkdir kp1_nude
cd kp1_nude
mkdir concat fastqc raw rsem trim

logout

#From sammas folder:
scp -r ss4cf$ /Volumes/JanesLab/RNA-seq\ raw\ data/08122019-Janes_Pool13-139986850* ss4cf@rivanna.hpc.virginia.edu:/scratch/ss4cf/kp1_nude/raw

#Navigate back to scratch folder:
ssh -Y ss4cf@rivanna.hpc.virginia.edu
#password

cd /scratch/ss4cf/kp1_nude
# Get the sample names
SAMPLIST=$(find ./raw/*/*/* | cut -f 5 -d "/" | cut -f 1 -d "_" | sort -u)

echo $SAMPLIST #to confirm the names

for SAMP in $SAMPLIST
do
R1=$(find ./raw/*/*/$SAMP*R1*.fastq.gz | sort)
R2=$(find ./raw/*/*/$SAMP*R2*.fastq.gz | sort)
cat $R1 > ./concat/${SAMP}_R1.fastq.gz
cat $R2 > ./concat/${SAMP}_R2.fastq.gz
done

#I usually go check that these files are actually in my concat folder

ls concat
#Two sets are named weird here, maybe a typo by Emily? fix in concat:
cd concat
mv KP1-Nude10c03_R1.fastq.gz KP1Nude-10c03_R1.fastq.gz
mv KP1-Nude10c03_R2.fastq.gz KP1Nude-10c03_R2.fastq.gz 
mv KP1-Nude-10c16_R1.fastq.gz KP1Nude-10c16_R1.fastq.gz
mv KP1-Nude-10c16_R2.fastq.gz KP1Nude-10c16_R2.fastq.gz 

cd..
# Trim adapter sequences - previously did this on the front end, but going to try to do it as a job
pwd
#/scratch/ss4cf/kp1_nude

#This is a dry run and produces a list of commands, test one to make sure there are no errors
find ./concat/*fastq.gz | cut -f 1 -d "_" | cut -f 3 -d "/" | sort | parallel -j 6 --dry-run "~/software/ea-utils/ea-utils-master/clipper/fastq-mcf -q 10 -t 0.01 -k 0 ~/adapters/nextera_adapters.fa ./concat/{}_R1.fastq.gz ./concat/{}_R2.fastq.gz -o ./trim/{}_R1.trim.fq.gz -o ./trim/{}_R2.trim.fq.gz > ./trim/{}.trim.log"
```

New Rivanna update issues: have to load gparallel module for above to work! also had to remake clipper/ea-utils as outlined in [Adapter trimming](Adapter trimming)


```bash
#make a trim.sh file:
nano trim.sh
--------
#!/bin/sh

module load gparallel

find ./concat/*fastq.gz | cut -f 1 -d "_" | cut -f 3 -d "/" | sort | parallel -j 6 "~/software/ea-utils/ea-utils-master/clipper/fastq-mcf -q 10 -t 0.01 -k 0 ~/adapters/nextera_adapters.fa ./concat/{}_R1.fastq.gz ./concat/{}_R2.fastq.gz -o ./trim/{}_R1.trim.fq.gz -o ./trim/{}_R2.trim.fq.gz > ./trim/{}.trim.log"

---------
sbatch -A janeslab -p standard -t 1:00:00 -n 1 --cpus-per-task 5 ./trim.sh
#can request up to 20, but then your job could wait pending resources etc...

cd trim
#Check number of reads before and after trimming:
SAMPLIST=$(find *log | grep -v _R | cut -f 1 -d "." | sort -u)
for SAMP in $SAMPLIST
do
RC=$(grep "Total reads" ${SAMP}.trim.log | cut -f 3 -d " ")
TR=$(grep "Too short" ${SAMP}.trim.log | cut -f 5 -d " ")
RL=$(( $RC - $TR ))
echo -e '|'$SAMP'|'$RC'|'$RL'|'
done

#Run FASTQC

cd ../fastqc
nano fastqcjob.sh
________________
#!/bin/bash

module load fastqc
module load gparallel

find ../trim/*fq.gz | parallel -j 6 "fastqc {}"
_____________

sbatch -A janeslab -p standard -t 1:00:00 -n 1 --cpus-per-task 10 ./fastqcjob.sh
cd ../trim
mv ./*fastqc* ../fastqc

```
#### Alignments with RSEM

```shell
cd ../rsem

nano submitjobs.sh
-------------
#!/bin/bash

SAMPLIST=$(find ../trim/*.fq.gz |cut -f 3 -d "/"|  cut -f 1 -d "_" | sort -u)

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

REF="/home/ss4cf/genomes/grcm38_cdna_82_ercc/grcm38ercc_rsem"

module load gparallel
module load gcc rsem

mkdir $SAMP

echo "Processing file $SAMP, starting at $(date +%F\ %T)"

rsem-calculate-expression -p 20 --bowtie2 --bowtie2-path /home/ss4cf/software/bowtie2/bowtie2-2.3.4.3 --single-cell-prior --paired-end ../trim/${SAMP}_R1.trim.fq.gz ../trim/${SAMP}_R2.trim.fq.gz $REF ./${SAMP}/$SAMP

echo "Finished processing file $SAMP at $(date +%F\ %T)"

--------------------

#Execute scripts:
chmod u+x *.sh #convert to executable
./submitjobs.sh #run
squeue -u ss4cf #check job queue/status
```
#### Multi QC on alignments

```shell
pwd
#/scratch/ss4cf/kp1_nude/rsem

# Load anaconda to use multiqc package which needs some python
module load anaconda
module load bioconda/py3.6
pip install --user multiqc

#Copy trim QC files to RSEM folder to consolidate
cp -r ../fastqc .
multiqc -f . #run multiQC

logout 

#Navigate to folder on your computer where you want to save it
cd /Users/ss4cf/Desktop/Lab/Dylan\ Schaff/RNAseq/KP1_nude 

#Copy multiQC results
scp -r ss4cf@rivanna.hpc.virginia.edu:/scratch/ss4cf/kp1_nude/rsem/multiqc* . 
#This includes a nice .html page to view results
```

#### Generating read count csv

```shell
logout
cd /Users/ss4cf/Desktop/Lab/Dylan\ Schaff/RNAseq/KP1_nude

scp ss4cf@rivanna.hpc.virginia.edu:/scratch/ss4cf/kp1_nude/rsem/*/*.isoforms.results .
```

#### HISAT2 alignments to human genome in parallel for comparison

Create environment:

```shell
export RNA_HOME=/scratch/ss4cf/kp1_nude
export RNA_DATA_DIR=/scratch/ss4cf/kp1_nude/raw
export RNA_DATA_TRIM_DIR=$RNA_HOME/trim
export RNA_REFS_DIR=/scratch/ss4cf/bcc10/hisat/refs
export RNA_REF_INDEX=$RNA_REFS_DIR/grch38
export RNA_REF_FASTA=$RNA_REF_INDEX.fa
export RNA_REF_GTF=$RNA_REFS_DIR/Homo_sapiens.GRCh38.86.gtf
export RNA_ALIGN_DIR=$RNA_HOME/hisat

cd hisat

nano submit_hisat.sh 
-----------
#!/bin/bash

SAMPLIST=$(find $RNA_DATA_TRIM_DIR/*.fq.gz |cut -f 6 -d "/" | cut -f 1 -d "_" | sort -u)

for SAMP in $SAMPLIST
do

sbatch -A janeslab -p standard -t 01:00:00 -n 1 --cpus-per-task 10 --wrap="sh ./hisat.sh $SAMP"
sleep 1 # wait 1 second between each job submission

done
-------------

nano hisat.sh 

----------------
#!/bin/bash

SAMP=$1

mkdir $SAMP

echo "Processing file $SAMP, starting at $(date +%F\ %T)"

/scratch/ss4cf/bcc10/hisat/student_tools/hisat2-2.1.0/hisat2 -p 20 --rg-id=${SAMP} -x $RNA_REF_INDEX/genome --dta -1 $RNA_DATA_TRIM_DIR/${SAMP}_R1.trim.fq.gz -2 $RNA_DATA_TRIM_DIR/${SAMP}_R2.trim.fq.gz -S ./${SAMP}/${SAMP}.sam

echo "Finished processing file $SAMP at $(date +%F\ %T)"

-----------------
chmod u+x *.sh #convert to executable


#run for all samples
./submit_hisat.sh
squeue -u ss4cf

#multiqc on these files as well
```
