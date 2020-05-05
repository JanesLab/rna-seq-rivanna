#### Make project directory

```bash
cd /scratch/mds5cg
mkdir opc_90dpi
cd opc_90dpi
mkdir raw
logout
```

#### Copy FASTQ sequencing files from server to Rivanna

```bash
scp -r "/Volumes/JanesLab/RNA-seq raw data/Janes-Pool9-115842727/." mds5cg@rivanna.hpc.virginia.edu:"/scratch/mds5cg/2019-01-opc_90dpi/raw/"
```

#### Concatenate files

```bash

ijob -A janeslab -p standard -c 2 -t 2:00:00

cd /scratch/mds5cg/opc_90dpi
mkdir concat

SAMPLIST=$(find ./raw/*/*/* | cut -f 5 -d "/" | cut -f 1 -d "_" | sort -u)
for SAMP in $SAMPLIST
do
R1=$(find ./raw/*/*/$SAMP*R1*.fastq.gz | sort)
R2=$(find ./raw/*/*/$SAMP*R2*.fastq.gz | sort)
cat $R1 > ./concat/${SAMP}_R1.fastq.gz
cat $R2 > ./concat/${SAMP}_R2.fastq.gz
done
```

#### Trim files

```bash
cd /scratch/mds5cg/opc_90dpi
mkdir trim

find ./concat/*fastq.gz | cut -f 1 -d "_" | cut -f 3 -d "/" | sort | parallel -j 24 "~/software/ea-utils/ea-utils-master/clipper/fastq-mcf -q 10 -t 0.01 -k 0 ~/adapters/nextera_adapters.fa ./concat/{}_R1.fastq.gz ./concat/{}_R2.fastq.gz -o ./trim/{}_R1.trim.fq.gz -o ./trim/{}_R2.trim.fq.gz > ./trim/{}.trim.log"
```

#### Run FASTQC

```bash
cd /scratch/mds5cg/opc_90dpi
mkdir fastqc

module load fastqc
#Run FastQC in parallel
find *fq.gz | parallel -j 24 "fastqc {}"
mv ./trim/*fastqc* ./fastqc

scancel -u mds5cg
```

#### Run RSEM alignments

```bash
cd /scratch/mds5cg/opc_90dpi
mkdir rsem
logout
```

#### Make the following files in text editor. Copy to rsem folder.

File1: **submitjobs.sh**

```bash
#!/bin/bash

SAMPLIST=$(find ../trim/*.fq.gz | cut -f 1 -d "_" | sort -u)

for SAMP in $SAMPLIST
do

sbatch -A janeslab -p standard -t 8:00:00 -n 1 --wrap="sh ./calculate_rsem.sh $SAMP"
sleep 1	# wait 1 second between each job submission
  
done
```

File2: **calculate_rsem.sh**

```bash
#!/bin/bash/

SAMP=$1

REF="/home/mds5cg/genomes/grcm38_cDNA_ercc/grcm38_cDNA_rsem"

module load rsem
module load bowtie2

mkdir $SAMP

echo "Processing file $SAMP, starting at $(date +%F\ %T)"

rsem-calculate-expression -p 24 --bowtie2 --single-cell-prior --paired-end ${SAMP}_R1.trim.fq.gz ${SAMP}_R2.trim.fq.gz $REF $SAMP

echo "Finished processing file $SAMP at $(date +%F\ %T)"

mv ../trim/${SAMP}.stat $SAMP
mv ../trim/${SAMP}.genes.results $SAMP
mv ../trim/${SAMP}*.bam $SAMP
```

#### Upload scripts

```bash
scp *.sh mds5cg@rivanna.hpc.virginia.edu:/scratch/mds5cg/opc_90dpi/rsem
```

#### Execute scripts

```bash
cd /scratch/mds5cg/opc_90dpi/rsem
chmod u+x *.sh
./submitjobs.sh
```

#### Move rsem files into sample folders and then into rsem folder

```bash
$SAMPLIST=$(find *.stat -type d | cut -f 1 -d ".")
for SAMP in $SAMPLIST; do mv ${SAMP}.stat $SAMP; done
for SAMP in $SAMPLIST; do mv ${SAMP}.genes.results $SAMP; done
for SAMP in $SAMPLIST; do mv ${SAMP}.isoforms.results $SAMP; done
for SAMP in $SAMPLIST; do mv ${SAMP}*.bam $SAMP; done
find * -type d -maxdepth 0 -exec mv {} ../rsem \;
logout
```

#### Download results

```bash
scp mds5cg@rivanna.hpc.virginia.edu:/scratch/mds5cg/opc_90dpi/rsem/*/*.isoforms.results .
```
