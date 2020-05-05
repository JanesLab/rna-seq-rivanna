### April 2019: Alignment of Dylan's second round of sequencing for KP1 liver mets in immune competent mice

#### Prepare genome (version 82 to match Alex's previous alignment)

```bash
## Download genome
mkdir grcm38_cdna_82
cd grcm38_cdna_82
wget ftp://ftp.ensembl.org/pub/release-82/fasta/mus_musculus/cdna/*
ls
# CHECKSUMS  Mus_musculus.GRCm38.cdna.abinitio.fa.gz  Mus_musculus.GRCm38.cdna.all.fa.gz  README

## Concatenate with ERCCs
gunzip Mus_musculus.GRCm38.cdna.all.fa.gz
cd ..
mkdir grcm38_cdna_82_ercc
cd grcm38_cdna_82_ercc
cat ../grcm38_cdna_82/Mus_musculus.GRCm38.cdna.all.fa ../ercc/ERCC92.fa > Mus_musculus.GRCm.38.cdna.all_plusERCC.fa
ls
#Mus_musculus.GRCm.38.cdna.all_plusERCC.fa

## Prepare with RSEM
pwd
#/home/ss4cf/genomes/grcm38_cdna_82_ercc
module load bowtie2
module load rsem
rsem-prepare-reference --bowtie2 --bowtie2-path /home/ss4cf/software/bowtie2/bowtie2-2.3.4.3 ./Mus_musculus.GRCm.38.cdna.all_plusERCC.fa grcm38ercc_rsem
ls
#grcm38ercc_rsem.1.bt2  grcm38ercc_rsem.4.bt2   grcm38ercc_rsem.n2g.idx.fa  grcm38ercc_rsem.seq             #Mus_musculus.GRCm.38.cdna.all_plusERCC.fa
#grcm38ercc_rsem.2.bt2  grcm38ercc_rsem.grp     grcm38ercc_rsem.rev.1.bt2   grcm38ercc_rsem.ti
#grcm38ercc_rsem.3.bt2  grcm38ercc_rsem.idx.fa  grcm38ercc_rsem.rev.2.bt2   grcm38ercc_rsem.transcripts.fa
```
#### Set up raw data 

``` bash
#Make a bunch of directories for later in **scratch** space
cd /scratch/ss4cf
mkdir kp1_impos 
cd kp1_impos
mkdir concat fastqc raw rsem trim
logout

#From sammas folder:
scp -r /Volumes/JanesLab/RNA-seq\ raw\ data/03042019-Janes__Pool10-126218093\ \(Dylan\'s\ KP1\ liver\ profiling\)/* ss4cf@rivanna.hpc.virginia.edu:/scratch/ss4cf/kp1_impos/raw

#Navigate back to scratch folder:
ssh -Y ss4cf@rivanna.hpc.virginia.edu
#password

cd /scratch/ss4cf/kp1_impos
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

cd concat
ls

# Trim adapter sequences - previously did this on the front end, but going to try to do it as a job
pwd
#/scratch/ss4cf/kp1_impos

#This is a dry run and produces a list of commands, test one to make sure there are no errors
find ./concat/*fastq.gz | cut -f 1 -d "_" | cut -f 3 -d "/" | sort | parallel -j 6 --dry-run "~/software/ea-utils/ea-utils-master/clipper/fastq-mcf -q 10 -t 0.01 -k 0 ~/adapters/nextera_adapters.fa ./concat/{}_R1.fastq.gz ./concat/{}_R2.fastq.gz -o ./trim/{}_R1.trim.fq.gz -o ./trim/{}_R2.trim.fq.gz > ./trim/{}.trim.log"


#make a trim.sh file:
nano trim.sh
--------
#!/bin/sh

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

find ../trim/ *fq.gz | parallel -j 6 "fastqc {}"
_____________

sbatch -A janeslab -p standard -t 1:00:00 -n 1 --cpus-per-task 10 ./fastqcjob.sh

cd ../trim
mv ./*fastqc* ../fastqc

```

#### Alignments with RSEM
```shell
cd ../rsem
cp ../../bcc5/rsem/*.sh . #previous job files written

#edit submitjobs.sh to request more cpus:
---------------
#!/bin/bash

SAMPLIST=$(find ../trim/*.fq.gz |cut -f 3 -d "/"|  cut -f 1 -d "_" | sort -u)

for SAMP in $SAMPLIST
do

sbatch -A janeslab -p standard -t 10:00:00 -n 1 --cpus-per-task 20 --wrap="sh ./calculate_rsem.sh $SAMP"
sleep 1 # wait 1 second between each job submission

done
-----------------

#Edit calculate_rsem.sh for reference genome

#!/bin/bash/

SAMP=$1

REF="/home/ss4cf/genomes/grcm38_cdna_82_ercc/grcm38ercc_rsem"

module load rsem

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

```bash
pwd
#/scratch/ss4cf/kp1_impos/rsem

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
scp -r ss4cf@rivanna.hpc.virginia.edu:/scratch/ss4cf/kp1_impos/rsem/multiqc* . 
#This includes a nice .html page to view results
```

#### Generating read count csv
```shell
logout
cd /Users/ss4cf/Desktop/Lab/Dylan\ Schaff/RNAseq/KP1_impos

scp ss4cf@rivanna.hpc.virginia.edu:/scratch/ss4cf/kp1_impos/rsem/*/*.isoforms.results .
```

#### Four samples had very poor alignment rates, try aligning these to the human genome: use the one from bcc5

```bash
cd /scratch/ss4cf/kp1_impos
mkdir rsem_h
cd rsem_h
cp ../../bcc5/rsem/*.sh . #previous job files written

nano submitjobs.sh
---------------
#!/bin/bash

SAMPLIST=$(find ../trim/KP1ImPos-{10c10,10c02,Pool19,Pool08}*.fq.gz |cut -f 3 -d "/"|  cut -f 1 -d "_" | sort -u)

for SAMP in $SAMPLIST
do

sbatch -A janeslab -p standard -t 10:00:00 -n 1 --wrap="sh ./calculate_rsem.sh $SAMP"
sleep 1	# wait 1 second between each job submission
  
done

-----------------
#Edit calculate_rsem.sh for reference genome

nano calculate_rsem.sh
-----------------
#!/bin/bash/

SAMP=$1

REF="/home/ss4cf/genomes/grch38_cdna_84_ercc/grch38_84_ercc_rsem"

module load rsem
module load bowtie2

mkdir $SAMP

echo "Processing file $SAMP, starting at $(date +%F\ %T)"

rsem-calculate-expression -p 20 --bowtie2 --single-cell-prior --paired-end ../trim/${SAMP}_R1.trim.fq.gz ../trim/${SAMP}_R2.trim.fq.gz $REF ./${SAMP}/$SAMP

echo "Finished processing file $SAMP at $(date +%F\ %T)"
---------------
chmod u+x *.sh #convert to executable
./submitjobs.sh #run
squeue -u ss4cf #check job queue/status

#Get a quick QC table of these alignments:

# Load anaconda to use multiqc package which needs some python
module load anaconda
module load bioconda/py3.6
pip install --user multiqc

#Copy trim QC files to RSEM folder to consolidate
multiqc -f . #run multiQC

logout 

#Navigate to folder on your computer where you want to save it
cd /Users/ss4cf/Desktop/Lab/Dylan\ Schaff/RNAseq/KP1_impos
mkdir human_align
cd human_align
 
#Copy multiQC results
scp -r ss4cf@rivanna.hpc.virginia.edu:/scratch/ss4cf/kp1_impos/rsem_h/multiqc* . 
#This includes a nice .html page to view results
```
##### alignment rates not improved by aligning to human transcriptome...should I try blasting random reads?

blasting random reads from head ../trim/KP1ImPos-10c02_R1.trim.fq.gz | gunzip yielded 100% matches to wheat....

#### HISAT2 alignment:
June 3rd 2019: trying against human genomic, in case there is contamination?

First this, then mouse genomic?

Create environment:

```bash
export RNA_HOME=/scratch/ss4cf/kp1_impos/hisat
export RNA_DATA_DIR=/scratch/ss4cf/kp1_impos/raw
cd hisat
mkdir trimmed
export RNA_DATA_TRIM_DIR=$RNA_HOME/trimmed
mkdir refs
export RNA_REFS_DIR=/scratch/ss4cf/bcc10/hisat/refs
cd $RNA_REFS_DIR
mkdir H_with_ERCC92
export RNA_REF_INDEX=$RNA_REFS_DIR/grch38
export RNA_REF_FASTA=$RNA_REF_INDEX.fa
export RNA_REF_GTF=$RNA_REFS_DIR/Homo_sapiens.GRCh38.86.gtf
export RNA_ALIGN_DIR=$RNA_HOME/alignments

```

Align already trimmed reads

```bash
cp ../trim/*.fq.gz ./trimmed

nano submit_hisat.sh 
-----------
#!/bin/bash

SAMPLIST=$(find $RNA_DATA_TRIM_DIR/*.fq.gz |cut -f 7 -d "/" | cut -f 1 -d "_" | sort -u)

for SAMP in $SAMPLIST
do

sbatch -A janeslab -p standard -t 10:00:00 -n 1 --cpus-per-task 10 --wrap="sh ./hisat.sh $SAMP"
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


nano testsubmit.sh

-----------
#!/bin/bash

SAMPLIST=$(find $RNA_DATA_TRIM_DIR/KP1ImPos-10c01*.fq.gz |cut -f 7 -d "/" | cut -f 1 -d "_" | sort -u)

for SAMP in $SAMPLIST
do

sbatch -A janeslab -p standard -t 10:00:00 -n 1 --cpus-per-task 10 --wrap="sh ./hisat.sh $SAMP"
sleep 1 # wait 1 second between each job submission

done
------------
chmod u+x *.sh #convert to executable
./testsubmit.sh #run


#run for all samples
./submit_hisat.sh
squeue -u ss4cf


#convert all .sam files to .bam files

nano submitsam.sh

-------

#!/bin/bash

SAMPLIST=$(find */*.sam |cut -f 2 -d "/"| cut -f 1 -d "." | sort -u)

for SAMP in $SAMPLIST
do

sbatch -A janeslab -p standard -t 10:00:00 -n 1 --cpus-per-task 10 --wrap="sh ./sam2bam.sh $SAMP"
sleep 1 # wait 1 second between each job submission

done

------

nano sam2bam.sh
----------------
#!/bin/bash

SAMP=$1

/scratch/ss4cf/bcc10/hisat/student_tools/samtools-1.9/samtools sort -@ 8 -o ./${SAMP}/${SAMP}.bam ./${SAMP}/${SAMP}.sam

-----------------

```

Two of these: 10c08 and Pool19 have >90% alignment to human genome...so there is definitely a contamination issue here.

The fact that these don't align to the human transcriptome via RSEM but do the genome makes me think it could be anything!? Like something fell into the reagents?

Do more QC to understand?

```bash
mkdir align_fastqc
fastqc ./*/*.bam
mv *fastqc.html align_fastqc/
mv *fastqc.zip align_fastqc/


```

Align KP1 samples to mouse genome via HISAT as well just for comparison/to see if it improves. Also, what is the level of overlap between intergenic regions of mouse and human genome?


```bash
export RNA_HOME=/scratch/ss4cf/kp1_impos/hisat
export RNA_DATA_DIR=/scratch/ss4cf/kp1_impos/raw
cd hisat
export RNA_DATA_TRIM_DIR=$RNA_HOME/trimmed
mkdir mouseref
export RNA_REFS_DIR=/scratch/ss4cf/kp1_impos/hisat/mouseref

#Download reference and matched index:
wget ftp://ftp.ensembl.org/pub/release-84/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz

wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grcm38.tar.gz

wget ftp://ftp.ensembl.org/pub/release-84/gtf/mus_musculus/Mus_musculus.GRCm38.84.gtf.gz 

#unzip untar all

export RNA_REF_INDEX=$RNA_REFS_DIR/grcm38
export RNA_REF_FASTA=$RNA_REF_INDEX.fa
export RNA_REF_GTF=$RNA_REFS_DIR/Mus_musculus.GRCm38.84.gtf
export RNA_ALIGN_DIR=$RNA_HOME/mousealignments
```
SUBMIT alignments:

```bash
nano submit_hisat.sh 
-----------
#!/bin/bash

SAMPLIST=$(find $RNA_DATA_TRIM_DIR/*.fq.gz |cut -f 7 -d "/" | cut -f 1 -d "_" | sort -u)

for SAMP in $SAMPLIST
do

sbatch -A janeslab -p standard -t 10:00:00 -n 1 --cpus-per-task 10 --wrap="sh ./hisat.sh $SAMP"
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


nano testsubmit.sh

-----------
#!/bin/bash

SAMPLIST=$(find $RNA_DATA_TRIM_DIR/KP1ImPos-10c01*.fq.gz |cut -f 7 -d "/" | cut -f 1 -d "_" | sort -u)

for SAMP in $SAMPLIST
do

sbatch -A janeslab -p standard -t 10:00:00 -n 1 --cpus-per-task 10 --wrap="sh ./hisat.sh $SAMP"
sleep 1 # wait 1 second between each job submission

done
------------
chmod u+x *.sh #convert to executable
./testsubmit.sh #run


#run for all samples
./submit_hisat.sh
squeue -u ss4cf




```

