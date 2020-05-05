#### Setting up Rivanna for sequencing 

I followed steps from `Sham BCC10 alignments.md` for the most part except where indicated.

#### Get sequencing data
I uploaded the relevant sequencing files from sammas value storage to a new scratch folder:

```bash
cd /scratch/ss4cf
mkdir bcc6 #we will combine both cartridges here
cd bcc6
#Make a bunch of directories for later (there is probably a better way to do this)
mkdir concat fastqc raw rsem trim

logout 
```
Now, I log in to the sammas storage on my computer as usual. Then, secure copy files from sammas to scratch folder.

This can take a while depending on number of samples (~10s/file).

```bash
#both runs
scp -r /Volumes/JanesLab/RNA-seq\ raw\ data/11132019-Janes_Pool14-144040901/* ss4cf@rivanna.hpc.virginia.edu:/scratch/ss4cf/bcc6/raw

```
#### Concatenate reads/sample
For each cartridge, we have paired end reads (R1 and R2) in 4 different lanes (L1-4). Next, we concatenate all the R1 and R2 reads for different lanes and both cartridges.

Log into Rivanna and navigate to your scratch folder for that account

```bash
pwd
#/scratch/ss4cf/bcc6 (path name)

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

ls concat

```
#### Trimming reads
Next, Illumina adapter sequences are trimmed from raw reads.

```bash
pwd
#/scratch/ss4cf/bcc6 (path name)

#This is a dry run and produces a list of commands, test one to make sure there are no errors
module load gparallel

find ./concat/*fastq.gz | cut -f 1 -d "_" | cut -f 3 -d "/" | sort | parallel -j 6 --dry-run "~/software/ea-utils/ea-utils-master/clipper/fastq-mcf -q 10 -t 0.01 -k 0 ~/adapters/nextera_adapters.fa ./concat/{}_R1.fastq.gz ./concat/{}_R2.fastq.gz -o ./trim/{}_R1.trim.fq.gz -o ./trim/{}_R2.trim.fq.gz > ./trim/{}.trim.log"


#make a trim.sh file:
nano trim.sh
--------
#!/bin/sh

module load gparallel

find ./concat/*fastq.gz | cut -f 1 -d "_" | cut -f 3 -d "/" | sort | parallel -j 6 "~/software/ea-utils/ea-utils-master/clipper/fastq-mcf -q 10 -t 0.01 -k 0 ~/adapters/nextera_adapters.fa ./concat/{}_R1.fastq.gz ./concat/{}_R2.fastq.gz -o ./trim/{}_R1.trim.fq.gz -o ./trim/{}_R2.trim.fq.gz > ./trim/{}.trim.log"

---------
sbatch -A janeslab -p standard -t 3:00:00 -n 1 --cpus-per-task 5 ./trim.sh
#can request up to 20, but then your job could wait pending resources etc...

```
#### QC of trimmed reads

```bash
#Check number of reads before and after trimming:

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

nano fastqcjob.sh
________________
#!/bin/bash

module load fastqc
module load gparallel

find trim/*fq.gz | parallel -j 6 "fastqc {}"
_____________

sbatch -A janeslab -p standard -t 1:00:00 -n 1 --cpus-per-task 10 ./fastqcjob.sh
cd ../trim
mv ./*fastqc* ../fastqc

```

#### Alignment with RSEM. 

This will take ~1h/million reads so adjust requested time per sample based on read-depth

```bash
cd ../rsem

File 1: submitjobs.sh

```bash
nano submitjobs.sh
-------------
#!/bin/bash

SAMPLIST=$(find ../trim/*.fq.gz |cut -f 3 -d "/"|  cut -f 1 -d "_" | sort -u)

for SAMP in $SAMPLIST
do

sbatch -A janeslab -p standard -t 10:00:00 -n 1 --cpus-per-task 5 --wrap="sh ./calculate_rsem.sh $SAMP"
sleep 1 # wait 1 second between each job submission

done
```

---------------
```

File 2: calculate_rsem.sh

This uses an updated version of bowtie2 that is different from the one on Rivanna (to match the Bioinformatics core). See [Using an updated version of Bowtie2](Using an updated version of Bowtie2). 

```bash

nano calculate_rsem.sh
------------
#!/bin/bash/

SAMP=$1

REF="/home/ss4cf/genomes/grch38_cdna_84_ercc/grch38_84_ercc_rsem"

module load rsem

module load gparallel
module load gcc rsem

mkdir $SAMP

echo "Processing file $SAMP, starting at $(date +%F\ %T)"

rsem-calculate-expression -p 20 --bowtie2 --bowtie2-path /home/ss4cf/software/bowtie2/bowtie2-2.3.4.3 --single-cell-prior --paired-end ../trim/${SAMP}_R1.trim.fq.gz ../trim/${SAMP}_R2.trim.fq.gz $REF ./${SAMP}/$SAMP

echo "Finished processing file $SAMP at $(date +%F\ %T)"
---------------
```

Execute scripts

**Note:** You can modify submitjobs.sh to only do 1 to test first

```bash
chmod u+x *.sh #convert to executable
./submitjobs.sh #run
squeue -u ss4cf #check job queue/status
```
#### Multi QC on alignments

```bash
# Load anaconda to use multiqc package which needs some python
module load anaconda
module load bioconda/py3.6
pip install --user multiqc

#Copy trim QC files to RSEM folder to consolidate
cp -r ../fastqc .
multiqc -f . #run multiQC

logout 

#Navigate to folder on your computer where you want to save it
cd /Users/ss4cf/Desktop/Lab/Data/BCC\ Samples/BCC_RNAseq/RSEM\ data/BCC6/BCC6_alignment

#Copy multiQC results
scp -r ss4cf@rivanna.hpc.virginia.edu:/scratch/ss4cf/bcc6/rsem/multiqc* . 
```

This includes a nice .html page to view results

#### Download results
Logout of Rivanna and navigate to the folder you want to save the resulting files 

```bash
#Save results on your personal computer:

logout
cd /Users⁩/ss4cf⁩/⁨Desktop⁩/⁨Lab⁩/⁨Data⁩/⁨BCC\ Samples⁩/⁨BCC_RNAseq⁩/⁨RSEM\ data⁩/⁨BCC6/BCC6_alignment/BCC6_alignment_results

scp ss4cf@rivanna.hpc.virginia.edu:/scratch/ss4cf/bcc6/rsem/*/*.isoforms.results .
```

To arrive at compiled list of all samples and genes refer to [CSV generation in R](CSV generation in R).

There are still some low alignment samples though not as severe as before


#### HISAT2 alignment: troubleshooting low alignment, looking for intronic coverage?

Create environment:
```bash
cd /scratch/ss4cf/bcc6/
mkdir hisat
export RNA_HOME=/scratch/ss4cf/bcc6/hisat
export RNA_DATA_DIR=/scratch/ss4cf/bcc6/raw
cd hisat
mkdir trimmed
mkdir alignments
export RNA_DATA_TRIM_DIR=/scratch/ss4cf/bcc6/trim
export RNA_REFS_DIR=/scratch/ss4cf/bcc10/hisat/refs #use same reference as BCC10
export RNA_REF_INDEX=$RNA_REFS_DIR/grch38
export RNA_REF_FASTA=$RNA_REF_INDEX.fa
export RNA_REF_GTF=$RNA_REFS_DIR/Homo_sapiens.GRCh38.86.gtf
export RNA_ALIGN_DIR=$RNA_HOME/alignments


#tools
PATH=/scratch/ss4cf/bcc10/hisat/student_tools/genePredToBed:/scratch/ss4cf/bcc10/hisat/student_tools/gtfToGenePred:/scratch/ss4cf/bcc10/hisat/student_tools/bedops_linux_x86_64-v2.4.35/bin:/scratch/ss4cf/bcc10/hisat/student_tools/samtools-1.9:/scratch/ss4cf/bcc10/hisat/student_tools/bam-readcount/bin:/scratch/ss4cf/bcc10/hisat/student_tools/hisat2-2.1.0:/scratch/ss4cf/bcc10/hisat/student_tools/stringtie-1.3.4d.Linux_x86_64:/scratch/ss4cf/bcc10/hisat/student_tools/gffcompare-0.10.6.Linux_x86_64:/scratch/ss4cf/bcc10/hisat/student_tools/htseq-release_0.11.0/scripts:/scratch/ss4cf/bcc10/hisat/student_tools/tophat-2.1.1.Linux_x86_64:/scratch/ss4cf/bcc10/hisat/student_tools/kallisto_linux-v0.44.0:/scratch/ss4cf/bcc10/hisat/student_tools/FastQC:/scratch/ss4cf/bcc10/hisat/student_tools/flexbar-3.4.0-linux:/scratch/ss4cf/bcc10/hisat/student_tools/regtools/build:/home/ubuntu/bin/bedtools2/bin:$PATH

export LD_LIBRARY_PATH=/scratch/ss4cf/bcc10/hisat/student_tools/flexbar-3.4.0-linux:$LD_LIBRARY_PATH

echo $PATH
```
Alignments:
I already adapter trimmed these for RSEM, so just use those trimmed reads

```bash
#in /scratch/ss4cf/bcc6/hisat/alignments

nano submit_hisat.sh 
-----------
#!/bin/bash

SAMPLIST=$(find $RNA_DATA_TRIM_DIR/*.fq.gz |cut -f 6 -d "/"|  cut -f 1 -d "_" | sort -u)

for SAMP in $SAMPLIST
do

sbatch -A janeslab -p standard -t 02:00:00 -n 1 --cpus-per-task 5 --wrap="sh ./hisat.sh $SAMP"
sleep 1 # wait 1 second between each job submission

done
------------

nano hisat.sh 

----------------
#!/bin/bash

SAMP=$1

mkdir $SAMP

echo "Processing file $SAMP, starting at $(date +%F\ %T)"

/scratch/ss4cf/bcc10/hisat/student_tools/hisat2-2.1.0/hisat2 -p 20 --rg-id=${SAMP} -x $RNA_REF_INDEX/genome --dta -1 $RNA_DATA_TRIM_DIR/${SAMP}_R1.trim.fq.gz -2 $RNA_DATA_TRIM_DIR/${SAMP}_R2.trim.fq.gz -S ./${SAMP}/${SAMP}.sam

echo "Finished processing file $SAMP at $(date +%F\ %T)"
-----------------

#modify submit to only be for 2 samples at first:

nano testsubmit.sh

-----------
#!/bin/bash

SAMPLIST=$(find $RNA_DATA_TRIM_DIR/BCC6-{Sample17,Control20}*.fq.gz |cut -f 6 -d "/"| cut -f 1 -d "_"| sort -u)

for SAMP in $SAMPLIST
do

sbatch -A janeslab -p standard -t 02:00:00 -n 1 --cpus-per-task 5 --wrap="sh ./hisat.sh $SAMP"
sleep 1 # wait 1 second between each job submission

done
------------
chmod u+x *.sh #convert to executable
./testsubmit.sh #run

#HMMM look into this more but Sample 12 now has 97% alignment...which indicates genomic/contamination...

#run for all samples
./submit_hisat.sh
squeue -u ss4cf

# try multiQC:

module load anaconda
module load bioconda/py3.6
pip install --user multiqc
cp -r /scratch/ss4cf/bcc6/fastqc .
multiqc -f . #run multiQC

#get sample names from slurm files

SAMPLIST=$(find slurm* | cut -f 1 -d ".")
for SAMP in $SAMPLIST
do
slurm=$(find ${SAMP}.out | cut -f 1 -d ".")
sample=$(grep "BCC6" ${SAMP}.out | cut -f 3 -d " " | cut -f 1 -d "," | head -1)
echo -e '|'$slurm'|'$sample'|'
done

#convert all .sam files to .bam files

nano submitsam.sh

-------

#!/bin/bash

SAMPLIST=$(find */*.sam |cut -f 2 -d "/"| cut -f 1 -d "." | sort -u)

for SAMP in $SAMPLIST
do

sbatch -A janeslab -p standard -t 2:00:00 -n 1 --cpus-per-task 5 --wrap="sh ./sam2bam.sh $SAMP"
sleep 1 # wait 1 second between each job submission

done

------

nano sam2bam.sh
----------------
#!/bin/bash

SAMP=$1

/scratch/ss4cf/bcc10/hisat/student_tools/samtools-1.9/samtools sort -@ 8 -o ./${SAMP}/${SAMP}.bam ./${SAMP}/${SAMP}.sam

-----------------
####Get HTseq raw counts and compare to RSEM output?

nano submit_htseq.sh
-------

#!/bin/bash

SAMPLIST=$(find ../*/*.sam |cut -f 2 -d "/"| cut -f 1 -d "." | sort -u)

for SAMP in $SAMPLIST
do

sbatch -A janeslab -p standard -t 2:00:00 -n 1 --cpus-per-task 5 --wrap="sh ./htseq_counts.sh $SAMP"
sleep 1 # wait 1 second between each job submission

done

------

nano htseq_counts.sh

----------------
#!/bin/bash

SAMP=$1

/scratch/ss4cf/bcc10/hisat/student_tools/htseq-release_0.11.0/scripts/htseq-count --format bam --order pos --mode intersection-strict --minaqual 1 --type exon --idattr gene_id $RNA_ALIGN_DIR/${SAMP}/${SAMP}.bam $RNA_REF_GTF >${SAMP}_gene.tsv

-----------------
