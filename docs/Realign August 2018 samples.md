
```bash
cd /scratch/mds5cg/example_all
SAMPLIST=$(find ./raw/*/*/* | cut -f 5 -d "/" | cut -f 1 -d "_" | sort -u)

for SAMP in $SAMPLIST
do
R1=$(find ./raw/*/*/$SAMP*R1*.fastq.gz | sort)
R2=$(find ./raw/*/*/$SAMP*R2*.fastq.gz | sort)
cat $R1 > ./concat/${SAMP}_R1.fastq.gz
cat $R2 > ./concat/${SAMP}_R2.fastq.gz
done

find ./concat/*fastq.gz | cut -f 1 -d "_" | cut -f 3 -d "/" | sort -u | parallel -j 24 --dry-run "~/software/ea-utils/ea-utils-master/clipper/fastq-mcf -q 10 -t 0.01 -k 0 ~/adapters/all-illumina-adapters-with-pfos-from-alison.fa ./concat/{}_R1.fastq.gz ./concat/{}_R2.fastq.gz -o ./trim_all/{}_R1.trim.fq.gz -o ./trim_all/{}_R2.trim.fq.gz > ./trim_all/{}.trim.log"

ijob -A janeslab -p standard -c 2 -t 2:00:00

find ./concat/*fastq.gz | cut -f 1 -d "_" | cut -f 3 -d "/" | sort -u | parallel -j 24 "~/software/ea-utils/ea-utils-master/clipper/fastq-mcf -q 10 -t 0.01 -k 0 ~/adapters/all-illumina-adapters-with-pfos-from-alison.fa ./concat/{}_R1.fastq.gz ./concat/{}_R2.fastq.gz -o ./trim_all/{}_R1.trim.fq.gz -o ./trim_all/{}_R2.trim.fq.gz > ./trim_all/{}.trim.log"

find ./trim_all/*fq.gz | parallel -j 24 "fastqc {}"

mv ./trim_all/*fastqc* ./fastqc
```

```bash
#!/bin/bash

REF="/home/mds5cg/genomes/grcm38_82/grcm38_82_ercc_bowtie2update/grcm38ercc_rsem"

SAMPLIST=$(find ../trim_alex/*bead1*fq.gz | cut -f 1-2 -d "_" | cut -f 3 -d "/" |  sort -u)

module load rsem
module load samtools
for SAMP in $SAMPLIST
do
rsem-calculate-expression -p 24 -bowtie2 --bowtie2-path /home/mds5cg/software/bowtie2/bowtie2-2.3.4.3 --single-cell-prior --paired-end ../trim_alex/${SAMP}_R1.trim.fq.gz ../trim_alex/${SAMP}_R2.trim.fq.gz $REF $SAMP
cp ../trim/${SAMP}.isoforms.results  .
done
```

First submission:

submitjobs.sh
```bash
#!/bin/bash

SAMPLIST=$(find ../trim_all/*.fq.gz | cut -f 1-2 -d "_" | sort -u)

for SAMP in $SAMPLIST
do

sbatch -A janeslab -p standard -t 8:00:00 -n 1 --wrap="sh ./calculate_rsem.sh $SAMP"
sleep 1 # wait 1 second between each job submission

done
```

calculate_rsem.sh
```bash
#!/bin/bash/

SAMP=$1

REF="/home/mds5cg/genomes/grcm38_82/grcm38_82_ercc_bowtie2update/grcm38ercc_rsem"

module load rsem

mkdir -p $SAMP

echo "Processing file $SAMP, starting at $(date +%F\ %T)"

rsem-calculate-expression -p 24 --bowtie2 --bowtie2-path /home/mds5cg/software/bowtie2/bowtie2-2.3.4.3 --single-cell-prior --paired-end ${SAMP}_R1.trim.fq.gz ${SAMP}_R2.trim.fq.gz $REF $SAMP

echo "Finished processing file $SAMP at $(date +%F\ %T)"
```

Five samples (KP-bead 2,4,6 and KP-nobead3,6) failed. Not enough time. Re-run with 12 hours.

submitjobs.sh
```bash
#!/bin/bash

SAMPLIST=(KP-bead2 KP-bead4 KP-bead6 KP-nobead3 KP-nobead6)

for SAMP in ${SAMPLIST[@]}
do

sbatch -A janeslab -p standard -t 12:00:00 -n 1 --wrap="sh ./calculate_rsem.sh ../trim_all/$SAMP"
sleep 1 # wait 1 second between each job submission

done
```

Next, I moved all files (*.bam *.isoforms.results *.genes.results *.stat) into the appropriately named sample folder, and moved those folders to the rsem folder.

```bash
# Download all .isoforms.results files
scp mds5cg@rivanna.hpc.virginia.edu:/scratch/mds5cg/example_all/rsem/*/*.isoforms.results .
```

Imported into R using the rsem2csv.R script and correlated with Alex's processing of the MADM-bead1 sample. Exact correlation (R^2=1).
