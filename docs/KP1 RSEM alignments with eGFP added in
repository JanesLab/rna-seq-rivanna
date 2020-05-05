Based on Sham's alignment of KP1 Nude mice Liver Mets August 2019 + extra cartridge December 2019 and other alignments just with an extra sequence added in.

Will have to add eGFP sequence and then prepare RSEM reference...


Got eGFP sequence from https://www.algosome.com/resources/common-sequences.html#gst
and created a fasta file for it "egfp.fa" in ss4cf/genomes

Then followed instructions on the tutorial to prepare the new concatenated reference.

Now go into KP1-Nude (try the latest 24 sample run for a quick comparison):

```shell
pwd
#/scratch/ss4cf/kp1_nude_seq2

mkdir rsem_egfp
cd rsem egfp


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

REF="/home/ss4cf/genomes/grcm38_cdna_82_ercc_egfp/grcm38_ercc_egfp_rsem"

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

Didn't really work, going to try with entire plx302-eGFP plasmid sequence.

