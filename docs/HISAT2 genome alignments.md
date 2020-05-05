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
mkdir H_with_ERCC92 #NEED TO CHANGE LATER, THERE ARE NO ERCC SEQUENCES CURRENTLY
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
