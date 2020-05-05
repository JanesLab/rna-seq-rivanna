Analyzing this similar to Ray's (Zong Lab) RNA-seq data 

#### Obtain raw data
Download FASTQ ftp links for RNA-seq data from https://www.ebi.ac.uk/ena/data/view/PRJEB30443

Copy text file containing these to scratch folder

```bash
pwd
#/scratch/ss4cf/ena_PRJEB30443/raw

#nano two files for doing this:

download.sh

SAMPLIST=$(cat wget_rna.txt | cut -f 1| tr -d '\r')

for SAMP in $SAMPLIST; do wget $SAMP; done

submitdownload.sh
#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH -A janeslab
#SBATCH -p standard
#SBATCH --nodes=2
#SBATCH --cpus-per-task=2


./download.sh

#then
chmod u+x *.sh
sbatch submitdownload.sh

```
### Parsing sample names
From Paper:
we performed RNA-sequencing on 22 WB1P tumors and 7 tumors from B1P mice injected with Lenti- Cre, and compared their expression profile to tumors from the KB1P mouse model and a mouse model of luminal breast cancer (WapCre;Cdh1F/F;PtenF/F, WEP; 25), using a three-gene signature that distinguishes the PAM50 subtypes

Kevin said we are interested in:

K14Cre;Brca1F/F;Trp53F/F (KB1P)

WapCre;Brca1F/F;Trp53F/F, WB1P

intraductal injection of lentiviral vectors22,23,24 expressing the Cre-recombinase (Lenti-Cre) in Brca1F/F;Trp53F/F (B1P)

Have to figure out which these are....DID IT!! had to use a combination of ENA information and SRA run table from NCBI 
https://www.ncbi.nlm.nih.gov/Traces/study/?WebEnv=NCID_1_122018220_130.14.22.76_5555_1565097239_2266731417_0MetA0_S_HStore&query_key=6

ONLY doing 22+10 repeats WB1P and 7 B1P for now: saved these runids in a txt file in raw directory - might move rest to a subdirectory

```bash
mkdir raw_all
mv *.gz raw_all
cd raw_all

SAMPLIST=$(cat ../runids.txt | cut -f 1| tr -d '\r')

for SAMP in $SAMPLIST; do file=$(find ./$SAMP*.gz); mv $file ../raw; 

done
#....this was VERY wrong...deleted everything and just re-wrote the last one

#redo...only download the 30 of interest, also submit as 30 different jobs

#this get just the run ids from the whole list
cat wget_rna.txt| find */*.gz |cut -f 2 -d "/"| cut -f 1 -d "." | sort -u

SAMPLIST=$(cat runids.txt | cut -f 1| tr -d '\r')

temp=$(grep "ERR3021148" wget_rna.txt)

#so i have to grep from this...?

#I made two files:

cat submitdownload.sh 
#!/bin/bash

SAMPLIST=$(cat runids.txt | cut -f 1| tr -d '\r')

for SAMP in $SAMPLIST
do
sbatch -A janeslab -p standard -t 1:00:00 -n 1 --cpus-per-task 2 --wrap="sh ./download.sh $SAMP"
sleep 1 # wait 1 second between each job submission
done

cat download.sh 
#!/bin/bash/

SAMP=$1
url=$(grep $SAMP wget_rna.txt |cut -f 1| tr -d '\r')
wget $url

#Trim reads
#dry run:
find ./raw/*fastq.gz | cut -f 2 -d "." | cut -f 3 -d "/" | sort | parallel -j 6 --dry-run "~/software/ea-utils/ea-utils-master/clipper/fastq-mcf -q 10 -t 0.01 -k 0 ~/adapters/nextera_adapters.fa ./raw/{}.fastq.gz -o ./trim/{}.trim.fq.gz > ./trim/{}.trim.log"

#make a trim.sh file:
nano trim.sh

--------
#!/bin/sh

find ./raw/*fastq.gz | cut -f 2 -d "." | cut -f 3 -d "/" | sort | parallel -j 6 "~/software/ea-utils/ea-utils-master/clipper/fastq-mcf -q 10 -t 0.01 -k 0 ~/adapters/nextera_adapters.fa ./raw/{}.fastq.gz -o ./trim/{}.trim.fq.gz > ./trim/{}.trim.log"

sbatch -A janeslab -p standard -t 1:00:00 -n 2 --cpus-per-task 5 ./trim.sh


#QC

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


#HISAT2 alignments: remember to modify for single end, check genome version used by Liz
export RNA_HOME=/scratch/ss4cf/ena_PRJEB30443
export RNA_DATA_DIR=/scratch/ss4cf/ena_PRJEB30443/raw
export RNA_DATA_TRIM_DIR=$RNA_HOME/trim
export RNA_REFS_DIR=/scratch/ss4cf/kp1_impos/hisat/mouseref #download for KP1
export RNA_REF_INDEX=$RNA_REFS_DIR/grcm38
export RNA_REF_FASTA=$RNA_REF_INDEX.fa
export RNA_REF_GTF=$RNA_REFS_DIR/Mus_musculus.GRCm38.84.gtf
export RNA_ALIGN_DIR=$RNA_HOME

#Submit HISAT

nano submit_hisat.sh 
-----------
#!/bin/bash

SAMPLIST=$(find $RNA_DATA_TRIM_DIR/*.fq.gz |cut -f 6 -d "/" | cut -f 1 -d "." | sort -u)

for SAMP in $SAMPLIST
do

sbatch -A janeslab -p standard -t 1:00:00 -n 1 --cpus-per-task 10 --wrap="sh ./hisat.sh $SAMP"
sleep 1 # wait 1 second between each job submission

done
-------------

nano hisat.sh 

----------------
#!/bin/bash

SAMP=$1

mkdir $SAMP

echo "Processing file $SAMP, starting at $(date +%F\ %T)"

/scratch/ss4cf/bcc10/hisat/student_tools/hisat2-2.1.0/hisat2 -p 20 --rg-id=${SAMP} -x $RNA_REF_INDEX/genome --dta -U $RNA_DATA_TRIM_DIR/${SAMP}.trim.fq.gz -S ./${SAMP}/${SAMP}.sam

echo "Finished processing file $SAMP at $(date +%F\ %T)"

-----------------
#input needs to be changed to -U for single/unpaired reads

nano testsubmit.sh

-----------
#!/bin/bash

SAMPLIST=$(find $RNA_DATA_TRIM_DIR/ERR3021220.trim.fq.gz |cut -f 6 -d "/" | cut -f 1 -d "." | sort -u)

for SAMP in $SAMPLIST
do

sbatch -A janeslab -p dev -t 01:00:00 -n 1 --cpus-per-task 5 --wrap="sh ./hisat.sh $SAMP"
sleep 1 # wait 1 second between each job submission

done
------------
#yay! this ran in ~10m

chmod u+x *.sh #convert to executable
./testsubmit.sh #run

#if this runs v fast, change time reqest

#run for all samples
./submit_hisat.sh
squeue -u ss4cf

#DONE! MultiQC:
module load anaconda
module load bioconda/py3.6
pip install --user multiqc
cp -r ../fastqc .
multiqc -f . #run multiQC

Convert all sam files to bam files (needed for StringTie):

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

```
#### Something I haven't done before: TPM estimation using StringTie which Liz did for Ray's data:

``` bash

#Attempting for ERR3021147 that finished 

#Stringtie (be in the stringie directory!!!) ...might even do this on the front end...tbd

SAMPLIST=$(find ../hisat/ERR*/*.bam|cut -f 3 -d "/"|cut -f 1 -d "."|sort -u)

for SAMP in $SAMPLIST
do
/scratch/ss4cf/bcc10/hisat/student_tools/stringtie-1.3.4d.Linux_x86_64/stringtie -G $RNA_REFS_DIR/Mus_musculus.GRCm38.95.gtf -p 20 -e -B -o ${SAMP}/transcripts.gtf -A ${SAMP}/gene_abundances.tsv $RNA_HOME/hisat/${SAMP}/${SAMP}.bam
done

#Testing for ERR3021147 that finished 
/scratch/ss4cf/bcc10/hisat/student_tools/stringtie-1.3.4d.Linux_x86_64/stringtie -G $RNA_REFS_DIR/Mus_musculus.GRCm38.95.gtf -p 20 -e -B -o ERR3021147/transcripts.gtf -A ERR3021147/gene_abundances.tsv $RNA_HOME/hisat/ERR3021147/ERR3021147.bam
#running....great! do below:

#Create expression matrix for all samples: 

#perl script from Griffith lab to combine files:
cd $RNA_HOME/hisat/stringtie/
wget https://raw.githubusercontent.com/griffithlab/rnabio.org/master/assets/scripts/stringtie_expression_matrix.pl
chmod +x stringtie_expression_matrix.pl

#get list of directorys
find ./ERR* -maxdepth 0 | cut -f 2 -d "/"|paste -sd, -

./stringtie_expression_matrix.pl --expression_metric=TPM --result_dirs='ERR3021147,ERR3021148,ERR3021149,ERR3021150,ERR3021151,ERR3021152,ERR3021153,ERR3021154,ERR3021155,ERR3021156,ERR3021157,ERR3021158,ERR3021159,ERR3021160,ERR3021161,ERR3021162,ERR3021163,ERR3021164,ERR3021165,ERR3021177,ERR3021178,ERR3021179,ERR3021208,ERR3021209,ERR3021210,ERR3021221,ERR3021222,ERR3021223,ERR3021224,ERR3021225,ERR3021226,ERR3021227,ERR3021228,ERR3021229,ERR3021230,ERR3021238,ERR3021239,ERR3021240' --transcript_matrix_file=transcript_tpm_all_samples.tsv --gene_matrix_file=gene_tpm_all_samples.tsv

#THIS WORKED ALSO, just modify for result_dirs!

```
