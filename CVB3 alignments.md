Log in to Rivanna (not shown). Make a folder for the CVB3 raw data. This is where we will upload files from the lab server.

```bash
mkdir -p /scratch/mds5cg/millie/0-raw
```

Remote log in to lab desktop (not shown), access the server via Command Prompt, and upload files to Rivanna

```shell
Y:
chdir "RNA-seq raw data\11022017-Janes_mRNASeq-51380329"
pscp -r * mds5cg@rivanna.hpc.virginia.edu:/scratch/mds5cg/millie/0-raw
```

I can disconnect from everything while this is happening. Uploading takes a few hours...

-----

I will now walk through the steps for the alignments, but note that all of these steps were completed by running one "master" script (provided at the end).

-----

Right now, we have multiple FASTQ files for each sample. This is because the samples were sequenced multiple times across 3 different cartridges to increase the read-depth (total number of reads per sample).

We want to combine ("concatenate") these files together so that we only have two files per sample: one for the "forward reads" and one for the "reverse reads", designated "R1" and "R2" in the file names.

```bash
# DO NOT RUN
cd /scratch/mds5cg/millie/

# Make a new directory for the concatenated files
mkdir -p /scratch/mds5cg/millie/1-cat

# Get the sample IDs from the file names, removing the directory information and the file extension
SAMPLES=$(find "0-raw" -name "*.fastq.gz" | rev | cut -f 1 -d "/" | rev | cut -f 1 -d "_" | sort -u)
# Output:
# 1C2-1
# 1C2-2
# 2E3-1
# 2E3-2
# ...

# Loop through the samples
for ISAMPLE in $SAMPLES
do
    # Gather all the R1 files and all the R2 files
    R1=$(find "./0-raw" -name "${ISAMPLE}*R1*.fastq.gz" | sort)
    R2=$(find "./0-raw" -name "${ISAMPLE}*R2*.fastq.gz" | sort)
    
    # Concatenate files and redirect the output to a new file
    # Simply calling "cat $R1" would print out the result rather than saving it to a file
    cat $R1 > "./1-cat/${ISAMPLE}_R1.fastq.gz"
    cat $R2 > "./1-cat/${ISAMPLE}_R2.fastq.gz"
done
```

Next, we want to check for any adapter sequences (used in the sequencing process) and "trim" them from the reads.

```bash
# DO NOT RUN
cd /scratch/mds5cg/millie/

# Make a new directory for the trimmed files (and the "log" files)
mkdir -p /scratch/mds5cg/millie/2-trim
mkdir -p /scratch/mds5cg/millie/2-trim/log

# Get the sample IDs from the file names, removing the directory information and the file extension
SAMPLES=$(find "0-raw" -name "*.fastq.gz" | rev | cut -f 1 -d "/" | rev | cut -f 1 -d "_" | sort -u)
# Output:
# 1C2-1
# 1C2-2
# 2E3-1
# 2E3-2
# ...

# Loop through the samples
for ISAMPLE in $SAMPLES
do
    ~/software/ea-utils/ea-utils-master/clipper/fastq-mcf -q 10 -t 0.01 -k 0 "/home/mds5cg/adapters/all-illumina-adapters-with-pfos-from-alison.fa" "./1-cat/${ISAMPLE}_R1.fastq.gz" "./1-cat/${ISAMPLE}_R2.fastq.gz" -o "./2-trim/${ISAMPLE}_R1.trim.fastq.gz" -o "./2-trim/${ISAMPLE}_R2.trim.fastq.gz" > "./2-trim/log/${ISAMPLE}.trim.log"
done
```

-----
Now, our files are ready for the HISAT2 alignments. But first, we need to obtain a reference genome. We will obtain a human genome reference from Ensembl, combine it with our SV40/CVB3 files, and build the HISAT2 index.

```bash
# DO NOT RUN
cd /scratch/mds5cg

# Make a folder for our reference genome
mkdir -p ./ref/human_sv40_cvb3
cd /scratch/mds5cg/ref/human_sv40_cvb3

# Download the FASTA file
wget ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
# Download the GTF annotation file
wget ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz

# ".gz" is a compressed file format. To eventually combine it with our custom sequences, we will need to expand it.
gunzip *.gz
```

Log out of Rivanna.

Upload custom sequences to Rivanna

[sv40_cvb3.fa](uploads/58cdbb20034362244025befde0447837/sv40_cvb3.fa)

[sv40_cvb3.gtf](uploads/fe4c3c62ab0d6cac7de3f451faa880c3/sv40_cvb3.gtf)

```bash
scp sv40_cvb3.fa mds5cg@rivanna.hpc.virginia.edu:/scratch/mds5cg/ref/human_sv40_cvb3
scp sv40_cvb3.gtf mds5cg@rivanna.hpc.virginia.edu:/scratch/mds5cg/ref/human_sv40_cvb3
```

Log in to Rivanna.

```bash
# DO NOT RUN
cd /scratch/mds5cg/ref/human_sv40_cvb3

# Concatenate files
cat Homo_sapiens.GRCh38.dna.primary_assembly.fa sv40_cvb3.fa > human_sv40_cvb3.fa
cat Homo_sapiens.GRCh38.99.gtf sv40_cvb3.gtf > human_sv40_cvb3.gtf

# Build HISAT2 index
module load gcc hisat2
hisat2-build human_sv40_cvb3.fa hisat2_human_sv40_cvb3
# This should create 8 new files
# hisat2_human_sv40_cvb3.1.ht2
# ...
# hisat2_human_sv40_cvb3.8.ht2
```

-----

Now we are ready to do the alignments.

```bash
# DO NOT RUN
cd /scratch/mds5cg/millie/

# Make a new directory for the output files
mkdir -p /scratch/mds5cg/millie/3-hisat2

# Get the sample IDs from the file names, removing the directory information and the file extension
SAMPLES=$(find "0-raw" -name "*.fastq.gz" | rev | cut -f 1 -d "/" | rev | cut -f 1 -d "_" | sort -u)
# Output:
# 1C2-1
# 1C2-2
# 2E3-1
# 2E3-2
# ...

# Create a variable "CORES", set equal to the number of processing cores we will use
CORES=1
# Create a variable "INDEX", set equal to the path of our HISAT2 index we previously built.
INDEX="/scratch/mds5cg/ref/human_sv40_cvb3/hisat2_human_sv40_cvb3"

module load gcc hisat2

# Loop through the samples
for ISAMPLE in $SAMPLES
do
    hisat2 -p $CORES --dta --rna-strandness RF --rg-id ${ISAMPLE} -x $INDEX -1 "./2-trim/${ISAMPLE}_R1.trim.fastq.gz" -2 "./2-trim/${ISAMPLE}_R2.trim.fastq.gz" -S "./3-hisat2/${ISAMPLE}.sam"
done
```

Explanation of options:

`-p $CORES` = number of processing cores to use in **p**arallel. Increasing will speed up the alignments

`--dta` = **d**ownstream **t**ranscriptome **a**ssembly. This is included because we're using Stringtie for assembly

`--rna-strandness RF` = because these are strand-specific reads from TruSeq strand-specific library

`--rg-id ${ISAMPLE}` = **r**ead **g**roup **ID**. It will print the sample ID name in the output header.

Now, we have to convert these output SAM files to BAM (binary, compressed) files, which will be fed into Stringtie.

```bash
# DO NOT RUN
cd /scratch/mds5cg/millie/

# Make a new directory for the output files
mkdir -p /scratch/mds5cg/millie/4-sam2bam

# Get the sample IDs from the file names, removing the directory information and the file extension
SAMPLES=$(find "0-raw" -name "*.fastq.gz" | rev | cut -f 1 -d "/" | rev | cut -f 1 -d "_" | sort -u)
# Output:
# 1C2-1
# 1C2-2
# 2E3-1
# 2E3-2
# ...

# Create a variable "CORES", set equal to the number of processing cores we will use
CORES=1

module load gcc samtools

# Loop through the samples
for ISAMPLE in $SAMPLES
do
    samtools sort -@ $CORES -o "./4-sam2bam/${SAMPLEID}.bam" "./3-hisat2/${SAMPLEID}.sam"
done
```

Finally, we will do the assembly using Stringtie.

```bash
# DO NOT RUN
cd /scratch/mds5cg/millie/

# Make a new directory for the output files
mkdir -p /scratch/mds5cg/millie/5-stringtie

# Get the sample IDs from the file names, removing the directory information and the file extension
SAMPLES=$(find "0-raw" -name "*.fastq.gz" | rev | cut -f 1 -d "/" | rev | cut -f 1 -d "_" | sort -u)
# Output:
# 1C2-1
# 1C2-2
# 2E3-1
# 2E3-2
# ...

# Create a variable "CORES", set equal to the number of processing cores we will use
CORES=1
# Create a variable "ANNOTATION", set equal to the path of our GTF annotation file.
INDEX="/scratch/mds5cg/ref/human_sv40_cvb3/human_sv40_cvb3.gtf"

module load gcc stringtie

# Loop through the samples
for ISAMPLE in $SAMPLES
do
    mkdir -p "./5-stringtie/${SAMPLEID}"
    stringtie -p $CORES -G $ANNOTATION -e -B -o "./5-stringtie/${SAMPLEID}/transcripts.gtf" -A "./5-stringtie/${SAMPLEID}/gene_abundances.tsv" "./4-sam2bam/${SAMPLEID}.bam"
done
```

Explanation of options:

`-p $CORES` = number of processing cores to use in **p**arallel. Increasing will speed up the alignments

`-G $ANNOTATION` = reference GTF annotation to use for assembly

`-e` = restrict the assembly to only known transcripts in the provided annotation

`-B` = save additional files for Ballgown (not necessary)

-----

The above steps can be combined into one master script, attached here.

[master.sh](uploads/b1e7b899199424dd63267adcc070bc16/master.sh)

-----

This creates one file per sample. However, it only includes TPM values. We want raw counts. To do this, we will use htseq.

```bash
# DO NOT RUN
cd /scratch/mds5cg/millie/

# Make a new directory for the output files
mkdir -p /scratch/mds5cg/millie/6-htseq

ANNOTATION="/scratch/mds5cg/ref/human_sv40_cvb3/human_sv40_cvb3.gtf"

SAMPLES=$(find ./4-sam2bam/*.bam | cut -f 3 -d "/" | cut -f 1 -d "." | sort -u)
for ISAMPLE in $SAMPLES
do
    htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id "./4-sam2bam/${ISAMPLE}.bam" $ANNOTATION > "./6-htseq/${ISAMPLE}.tsv"
done
```

We next want to join these files together into one big table where rows = genes and columns = samples.

First, we need this utility function (sourced from Stack overflow: https://stackoverflow.com/a/17649149) that will allow us to join more than 2 files at a time:

File name: `multijoin.sh`

```bash
#!/bin/sh

# multijoin - join multiple files

join_rec() {
    if [ $# -eq 1 ]; then
        join - "$1"
    else
        f=$1; shift
        join - "$f" | join_rec "$@"
    fi
}

if [ $# -le 2 ]; then
    join "$@"
else
    f1=$1; f2=$2; shift 2
    join "$f1" "$f2" | join_rec "$@"
fi
```

Now, combine the files using that utility function.

```bash
# DO NOT RUN
cd /scratch/mds5cg/millie/

SAMPLES=$(find ./6-htseq/*)

./multijoin.sh $SAMPLES > counts.tsv
# Output:
# CVB3_antisense   1     0     0     0     ...
# CVB3_sense       0     0     0     0     ...
# ENSG00000000003  978   1174  651   713   ...
# ENSG00000000005  0     0     0     0     ...
# ...              ...   ...   ...   ...   ...
```

Note that this file is missing a header. We want to include one, where the column titles are the sample IDs.

```bash
SAMPLELIST=$(find "7-htseq" -name "*.tsv" | cut -f 2 -d "/" | cut -f 1 -d "." | sort -u | paste -sd"\t")
echo -e "GeneID\t$SAMPLELIST" > header.tsv

# Replace the spaces with tab characters in the counts file
tr " " "\t" < counts.tsv > counts_mod.tsv

cat header.tsv counts_mod.tsv > counts_table.tsv
```
