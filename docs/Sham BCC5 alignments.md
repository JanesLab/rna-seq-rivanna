#### Setting up Rivanna for sequencing 

I followed instructions from the tutorial and downloaded the adapter file.

**Note:** When downloading a reference transcriptome, make sure to change grc**m** to grc**h** depending on whether you're aligning mouse or human samples.

**Note:** If comparing to previous runs aligned by the Bioinformatics core, ensure that you're using the same version of the reference genome and bowtie (aligner).


#### Get sequencing data
I uploaded the relevant sequencing files from sammas value storage to a new scratch folder:

```bash
cd /scratch/ss4cf
mkdir bcc5both #I'm combining two different runs to compare to a single run
cd bcc5both
#Make a bunch of directories for later (there is probably a better way to do this)
mkdir concat
mkdir fastqc
mkdir raw
mkdir rsem
mkdir trim

logout 
```
Now, I log in to the sammas storage on my computer as usual. Then, secure copy files from sammas to scratch folder. For this project, the two cartridges are in separate folders because they were run on different days. 

This can take a while depending on number of samples (~10s/file).

```bash
#run 1
scp -r /Volumes/JanesLab/RNA-seq\ raw\ data/09282018-Janes_Pool7-97687590\ \(Sham\'s\ BCC5\ profiling\)/* ss4cf@rivanna.hpc.virginia.edu:/scratch/ss4cf/bcc5both/raw

#run 2
scp -r /Volumes/JanesLab/RNA-seq\ raw\ data/09282018-Janes_Pool7_01-118819701\ \(Sham\'s\ BCC5\ profiling\,\ 2nd\ cartridge\)/* ss4cf@rivanna.hpc.virginia.edu:/scratch/ss4cf/bcc5both/raw
```

#### Concatenate reads/sample
For each cartridge, we have paired end reads (R1 and R2) in 4 different lanes (L1-4). Next, we concatenate all the R1 and R2 reads for different lanes and both cartridges.

Log into Rivanna and navigate to your scratch folder for that account
```bash
pwd
#/scratch/ss4cf/bcc5both (path name)

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
#/scratch/ss4cf/bcc5both (path name)

#This is a dry run and produces a list of commands, test one to make sure there are no errors
find ./concat/*fastq.gz | cut -f 1 -d "_" | cut -f 3 -d "/" | sort | parallel -j 6 --dry-run "~/software/ea-utils/ea-utils-master/clipper/fastq-mcf -q 10 -t 0.01 -k 0 ~/adapters/nextera_adapters.fa ./concat/{}_R1.fastq.gz ./concat/{}_R2.fastq.gz -o ./trim/{}_R1.trim.fq.gz -o ./trim/{}_R2.trim.fq.gz > ./trim/{}.trim.log"


find ./concat/*fastq.gz | cut -f 1 -d "_" | cut -f 3 -d "/" | sort | parallel -j 6 "~/software/ea-utils/ea-utils-master/clipper/fastq-mcf -q 10 -t 0.01 -k 0 ~/adapters/nextera_adapters.fa ./concat/{}_R1.fastq.gz ./concat/{}_R2.fastq.gz -o ./trim/{}_R1.trim.fq.gz -o ./trim/{}_R2.trim.fq.gz > ./trim/{}.trim.log"

#zipped files are necessary for bowtie
```
#### QC of trimmed reads

```bash
#Check number of reads before and after trimming:

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
cp ../../bcc5/rsem/*.sh . # I had written job files for a previous alignment, so I'm moving them into this current folder
```
We will need two executable files

**Note:** These files can be written in text edit on your local computer and then copied (scp) to Rivanna and then edited using vim/touch, see: [Writing scripts in Rivanna](Writing scripts in Rivanna).

File 1: submitjobs.sh

```bash
#!/bin/bash

SAMPLIST=$(find ../trim/*.fq.gz |cut -f 3 -d "/"|  cut -f 1 -d "_" | sort -u)

for SAMP in $SAMPLIST
do

sbatch -A janeslab -p standard -t 10:00:00 -n 1 --wrap="sh ./calculate_rsem.sh $SAMP"
sleep 1 # wait 1 second between each job submission

done
```

File 2: calculate_rsem.sh

This uses an updated version of bowtie2 that is different from the one on Rivanna (to match the Bioinformatics core). See [Using an updated version of Bowtie2](Using an updated version of Bowtie2). 

```bash
#!/bin/bash/

SAMP=$1

REF="/home/ss4cf/genomes/grch38_cdna_84_ercc/grch38_84_ercc_rsem"

module load rsem

mkdir $SAMP

echo "Processing file $SAMP, starting at $(date +%F\ %T)"

rsem-calculate-expression -p 24 --bowtie2 --bowtie2-path /home/ss4cf/software/bowtie2/bowtie2-2.3.4.3 --single-cell-prior --paired-end ../trim/${SAMP}_R1.trim.fq.gz ../trim/${SAMP}_R2.trim.fq.gz $REF ./${SAMP}/$SAMP

echo "Finished processing file $SAMP at $(date +%F\ %T)"
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
cd /Users/ss4cf/Desktop/Lab/Data/BCC\ Samples/BCC_RNAseq/RSEM\ data/BCC5/BCC5_alignment

#Copy multiQC results
scp -r ss4cf@rivanna.hpc.virginia.edu:/scratch/ss4cf/bcc5both/rsem/multiqc* . 
```

This includes a nice .html page to view results

#### Download results
Logout of Rivanna and navigate to the folder you want to save the resulting files 

```bash
#Save results on your personal computer:

logout
cd /Users⁩/ss4cf⁩/⁨Desktop⁩/⁨Lab⁩/⁨Data⁩/⁨BCC\ Samples⁩/⁨BCC_RNAseq⁩/⁨RSEM\ data⁩/⁨BCC5/BCC5_alignment/BCC5_alignment_results

scp ss4cf@rivanna.hpc.virginia.edu:/scratch/ss4cf/bcc5both/rsem/*/*.isoforms.results .
```

To arrive at compiled list of all samples and genes refer to [CSV generation in R](CSV generation in R).
