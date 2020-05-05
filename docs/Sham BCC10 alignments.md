#### Setting up Rivanna for sequencing 

I followed instructions from the tutorial and downloaded the adapter file.

**Note:** When downloading a reference transcriptome, make sure to change grc**m** to grc**h** depending on whether you're aligning mouse or human samples.

**Note:** If comparing to previous runs aligned by the Bioinformatics core, ensure that you're using the same version of the reference genome and bowtie (aligner).


#### Get sequencing data
I uploaded the relevant sequencing files from sammas value storage to a new scratch folder:

```bash
cd /scratch/ss4cf
mkdir bcc10 #we will combine both cartridges here
cd bcc10
#Make a bunch of directories for later (there is probably a better way to do this)
mkdir concat
mkdir fastqc
mkdir raw
mkdir rsem
mkdir trim

logout 
```
Now, I log in to the sammas storage on my computer as usual. Then, secure copy files from sammas to scratch folder.

This can take a while depending on number of samples (~10s/file).

```bash
#both runs
scp -r /Volumes/JanesLab/RNA-seq\ raw\ data/05212019-Janes_Pool11-133463330/* ss4cf@rivanna.hpc.virginia.edu:/scratch/ss4cf/bcc10/raw

```
#### Concatenate reads/sample
For each cartridge, we have paired end reads (R1 and R2) in 4 different lanes (L1-4). Next, we concatenate all the R1 and R2 reads for different lanes and both cartridges.

Log into Rivanna and navigate to your scratch folder for that account

```bash
pwd
#/scratch/ss4cf/bcc10 (path name)

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

```
#### Trimming reads
Next, Illumina adapter sequences are trimmed from raw reads.

```bash
pwd
#/scratch/ss4cf/bcc10 (path name)

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
module load fastqc
#Run FastQC in parallel
find *fq.gz | parallel -j 6 "fastqc {}"
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

sbatch -A janeslab -p standard -t 10:00:00 -n 1 --cpus-per-task 20 --wrap="sh ./calculate_rsem.sh $SAMP"
sleep 1 # wait 1 second between each job submission

done

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
cd /Users/ss4cf/Desktop/Lab/Data/BCC\ Samples/BCC_RNAseq/RSEM\ data/BCC10/BCC10_alignment

#Copy multiQC results
scp -r ss4cf@rivanna.hpc.virginia.edu:/scratch/ss4cf/bcc10/rsem/multiqc* . 
```

This includes a nice .html page to view results

#### Download results
Logout of Rivanna and navigate to the folder you want to save the resulting files 

```bash
#Save results on your personal computer:

logout
cd /Users⁩/ss4cf⁩/⁨Desktop⁩/⁨Lab⁩/⁨Data⁩/⁨BCC\ Samples⁩/⁨BCC_RNAseq⁩/⁨RSEM\ data⁩/⁨BCC10/BCC10_alignment/BCC10_alignment_results

scp ss4cf@rivanna.hpc.virginia.edu:/scratch/ss4cf/bcc10/rsem/*/*.isoforms.results .
```

To arrive at compiled list of all samples and genes refer to [CSV generation in R](CSV generation in R).

#### Troubleshooting low alignment rates

#pretty sure this is an actual sample issue, but going to maybe try a couple of alignment troubleshooting steps to see if i can figure this out

```bash
#Look at random sequences:
head BCC10-Sample12_R1.trim.fq.gz | gunzip
#blast these
```
Tried a couple of these, and it seems like again synthetic chromosomes etc - wondering if its genomic/DNA contamination after all?

Try trimming differently, more aggressively to get rid of AL1 primers?

Do this for BCC10-Sample12 which had the lowest alignment rate?

```bash
pwd
#/scratch/ss4cf/bcc10
mkdir retrim
#unzip trimmed files
gunzip trim/BCC10-Sample12_R1.trim.fq.gz
gunzip trim/BCC10-Sample12_R2.trim.fq.gz

#clip 15bp from ends - based on stuff Alex tried in the past. If this helps, can investigate further

#download seqtk
git clone https://github.com/lh3/seqtk.git
cd seqtk
make
export seqtk=/scratch/ss4cf/bcc10/retrim/seqtk/seqtk

#clip
$seqtk trimfq -b 15 ../trim/BCC10-Sample12_R1.trim.fq > BCC10-Sample12_R1.trim.clip.fq

$seqtk trimfq -b 15 ../trim/BCC10-Sample12_R2.trim.fq > BCC10-Sample12_R2.trim.clip.fq

#zip
gzip *.fq
```

Re-align these trimmed+clipped reads:

File 1: submitjobs.sh

```bash
nano submitjobs.sh
-------------
#!/bin/bash

SAMPLIST=$(find *.fq.gz |cut -f 3 -d "/"|  cut -f 1 -d "_" | sort -u)

for SAMP in $SAMPLIST
do

sbatch -A janeslab -p standard -t 4:00:00 -n 1 --cpus-per-task 10 --wrap="sh ./calculate_rsem.sh $SAMP"
sleep 1 # wait 1 second between each job submission

done

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

mkdir $SAMP

echo "Processing file $SAMP, starting at $(date +%F\ %T)"

rsem-calculate-expression -p 20 --bowtie2 --bowtie2-path /home/ss4cf/software/bowtie2/bowtie2-2.3.4.3 --single-cell-prior --paired-end ./${SAMP}_R1.trim.clip.fq.gz ./${SAMP}_R2.trim.clip.fq.gz $REF ./$SAMP

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
##### Re-trim front AND back and try once more
```bash
pwd
#/scratch/ss4cf/bcc10
mkdir retrim
#unzip trimmed files
gunzip trim/BCC10-Sample12_R1.trim.fq.gz
gunzip trim/BCC10-Sample12_R2.trim.fq.gz

#clip 10bp from both ends - https://github.com/lh3/seqtk

#download seqtk
git clone https://github.com/lh3/seqtk.git
cd seqtk
make
export seqtk=/scratch/ss4cf/bcc10/retrim/seqtk/seqtk

#clip
$seqtk trimfq -b 10 -e 10 ../trim/BCC10-Sample12_R1.trim.fq > BCC10-Sample12_R1.trim.clip.fq

$seqtk trimfq -b 10 -e 10 ../trim/BCC10-Sample12_R2.trim.fq > BCC10-Sample12_R2.trim.clip.fq

#zip
gzip *.fq

chmod u+x *.sh #convert to executable
./submitjobs.sh #run
squeue -u ss4cf #check job queue/status

```
#### Try aligning 1 good and 1 bad sample R1 reads only
```bash
mkdir singleend
cd singleend

gzip ../trim/*.fq

nano submitjobs.sh
-------------
#!/bin/bash

SAMPLIST=$(find ../trim/BCC10-{Sample12,Sample02}*.fq.gz |cut -f 3 -d "/"| cut -f 1 -d "_" | sort -u)

for SAMP in $SAMPLIST
do

sbatch -A janeslab -p standard -t 4:00:00 -n 1 --cpus-per-task 10 --wrap="sh ./calculate_rsem.sh $SAMP"

sleep 1 # wait 1 second between each job submission

done

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

mkdir $SAMP

echo "Processing file $SAMP, starting at $(date +%F\ %T)"

rsem-calculate-expression -p 20 --bowtie2 --bowtie2-path /home/ss4cf/software/bowtie2/bowtie2-2.3.4.3 --single-cell-prior ../trim/${SAMP}_R1.trim.fq.gz $REF ./${SAMP}/$SAMP

echo "Finished processing file $SAMP at $(date +%F\ %T)"
---------------

chmod u+x *.sh #convert to executable
./submitjobs.sh #run

```

#### Confirm that R1 and R2 reads are paired correctly?
```bash
zcat BCC10-Sample12_R2.trim.clip.fq.gz | grep "@NS500146:321:HKGWTAFXY:1:11101:19443:1083"
#@NS500146:321:HKGWTAFXY:1:11101:19443:1083 2:N:0:GGAGCTAC+NATGCAGT


zcat BCC10-Sample12_R1.trim.clip.fq.gz | grep "@NS500146:321:HKGWTAFXY:1:11101:19443:1083"
#@NS500146:321:HKGWTAFXY:1:11101:19443:1083 1:N:0:GGAGCTAC+NATGCAGT

#these seem identical, but maybe worth just repeating the concatenation? the reads are sorted though so this shouldn't be an issue.

```
This is the same...what does it mean?

#### HISAT2 alignment:

Might be worth trying a non-RSEM strategy for alignment?
reference free or HISAT?


Create environment:
```bash
export RNA_HOME=/scratch/ss4cf/bcc10/hisat
export RNA_DATA_DIR=/scratch/ss4cf/bcc10/raw
cd hisat
mkdir trimmed
export RNA_DATA_TRIM_DIR=$RNA_HOME/trimmed
mkdir refs
export RNA_REFS_DIR=$RNA_HOME/refs
cd $RNA_REFS_DIR
mkdir H_with_ERCC92
export RNA_REF_INDEX=$RNA_REFS_DIR/grch38
export RNA_REF_FASTA=$RNA_REF_INDEX.fa
export RNA_REF_GTF=$RNA_REFS_DIR/Homo_sapiens.GRCh38.86.gtf
export RNA_ALIGN_DIR=$RNA_HOME/alignments

```

#Download all tools:https://rnabio.org/module-00-setup/0000/08/01/Installation/

#Issue with the HTseq-count installment right now because of double python copies i think...
#Home directory was too full and thats where the python files themselves were...so moved some stuff to stratch and it worked!

#Also had issues with maybe flexbar and RSeQC...will revisit...

```bash
PATH=$RNA_HOME/student_tools/genePredToBed:$RNA_HOME/student_tools/gtfToGenePred:$RNA_HOME/student_tools/bedops_linux_x86_64-v2.4.35/bin:$RNA_HOME/student_tools/samtools-1.9:$RNA_HOME/student_tools/bam-readcount/bin:$RNA_HOME/student_tools/hisat2-2.1.0:$RNA_HOME/student_tools/stringtie-1.3.4d.Linux_x86_64:$RNA_HOME/student_tools/gffcompare-0.10.6.Linux_x86_64:$RNA_HOME/student_tools/htseq-release_0.11.0/scripts:$RNA_HOME/student_tools/tophat-2.1.1.Linux_x86_64:$RNA_HOME/student_tools/kallisto_linux-v0.44.0:$RNA_HOME/student_tools/FastQC:$RNA_HOME/student_tools/flexbar-3.4.0-linux:$RNA_HOME/student_tools/regtools/build:/home/ubuntu/bin/bedtools2/bin:$PATH

export LD_LIBRARY_PATH=$RNA_HOME/student_tools/flexbar-3.4.0-linux:$LD_LIBRARY_PATH

echo $PATH
```

Get a genome reference:

```bash
cd $RNA_HOME
echo $RNA_REFS_DIR
mkdir -p $RNA_REFS_DIR
wget ftp://ftp.ensembl.org/pub/release-86/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz

#annotations
echo $RNA_REFS_DIR
cd $RNA_REFS_DIR
wget ftp://ftp.ensembl.org/pub/release-86/gtf/homo_sapiens/Homo_sapiens.GRCh38.86.gtf.gz
echo $RNA_REF_GTF

gunzip *.gz

#Download a pre-indexed genome:
cd $RNA_REFS_DIR
wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch38.tar.gz
tar -zxvf grch38.tar.gz
#unzips into a bunch of files
```

Alignments:
I already adapter trimmed these for RSEM, so just use those trimmed reads

```bash
#in /scratch/ss4cf/bcc10/hisat/alignments

nano submit_hisat.sh 
-----------
#!/bin/bash

SAMPLIST=$(find $RNA_DATA_TRIM_DIR/*.fq.gz |cut -f 1 -d "_"| cut -f 7 -d "/" | sort -u)

for SAMP in $SAMPLIST
do

sbatch -A janeslab -p standard -t 10:00:00 -n 1 --cpus-per-task 10 --wrap="sh ./hisat.sh $SAMP"
sleep 1 # wait 1 second between each job submission

done
------------

nano hisat.sh 

----------------
#!/bin/bash

SAMP=$1

mkdir $SAMP

echo "Processing file $SAMP, starting at $(date +%F\ %T)"

$RNA_HOME/student_tools/hisat2-2.1.0/hisat2 -p 20 --rg-id=${SAMP} -x $RNA_REF_INDEX/genome --dta -1 $RNA_DATA_TRIM_DIR/${SAMP}_R1.trim.fq.gz -2 $RNA_DATA_TRIM_DIR/${SAMP}_R2.trim.fq.gz -S ./${SAMP}/${SAMP}.sam

echo "Finished processing file $SAMP at $(date +%F\ %T)"
-----------------

#modify submit to only be for 2 samples at first:

nano testsubmit.sh

-----------
#!/bin/bash

SAMPLIST=$(find $RNA_DATA_TRIM_DIR/BCC10-{Sample12,Sample02}*.fq.gz |cut -f 1 -d "_"| cut -f 7 -d "/" | sort -u)

for SAMP in $SAMPLIST
do

sbatch -A janeslab -p standard -t 10:00:00 -n 1 --cpus-per-task 10 --wrap="sh ./hisat.sh $SAMP"
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
cp -r /scratch/ss4cf/bcc10/fastqc .
multiqc -f . #run multiQC


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

$RNA_HOME/student_tools/samtools-1.9/samtools sort -@ 8 -o ./${SAMP}/${SAMP}.bam ./${SAMP}/${SAMP}.sam

-----------------

```

#Get HTseq raw counts and compare to RSEM output?

```bash

nano submit_htseq.sh
-------

#!/bin/bash

SAMPLIST=$(find ../*/*.sam |cut -f 2 -d "/"| cut -f 1 -d "." | sort -u)

for SAMP in $SAMPLIST
do

sbatch -A janeslab -p standard -t 10:00:00 -n 1 --cpus-per-task 10 --wrap="sh ./htseq_counts.sh $SAMP"
sleep 1 # wait 1 second between each job submission

done

------

nano htseq_counts.sh

----------------
#!/bin/bash

SAMP=$1

$RNA_HOME/student_tools/htseq-release_0.11.0/scripts/htseq-count --format bam --order pos --mode intersection-strict --minaqual 1 --type exon --idattr gene_id $RNA_ALIGN_DIR/${SAMP}/${SAMP}.bam $RNA_REF_GTF >${SAMP}_gene.tsv

-----------------




```
