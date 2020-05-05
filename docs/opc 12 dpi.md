Make scratch directory
```bash
cd /scratch/mds5cg
mkdir opc12
cd opc12
mkdir raw
logout
```

Copy files over from first cartridge
```bash
scp -r "/Volumes/JanesLab/RNA-seq raw data/10162019-Janes_Matt-opc12-142551423/FASTQ_Generation_2019-10-16_10_09_20Z-194261099/." mds5cg@rivanna.hpc.virginia.edu:/scratch/mds5cg/opc12/raw/
```

Concatenate
```bash
ijob -A janeslab -p standard -c 2 -t 2:00:00

SAMPLIST=$(find ./0-raw/*/* | cut -f 4 -d "/" | cut -f 1 -d "_" | sort -u)
for SAMP in $SAMPLIST
do
R1=$(find ./0-raw/*/$SAMP*R1*.fastq.gz | sort)
R2=$(find ./0-raw/*/$SAMP*R2*.fastq.gz | sort)
cat $R1 > ./1-concat/${SAMP}_R1.fastq.gz
cat $R2 > ./1-concat/${SAMP}_R2.fastq.gz
done
```

Trim
```bash
module load gcc gparallel
find ./1-concat/*fastq.gz | cut -f 3 -d "/" | cut -f 1 -d "_" | sort | parallel -j 24 "~/software/ea-utils/ea-utils-master/clipper/fastq-mcf -q 10 -t 0.01 -k 0 ~/adapters/nextera_adapters.fa ./1-concat/{}_R1.fastq.gz ./1-concat/{}_R2.fastq.gz -o ./2-trim/{}_R1.trim.fq.gz -o ./2-trim/{}_R2.trim.fq.gz > ./2-trim/{}.trim.log"
```

Align
```
#!/bin/bash
# submit_array.sh

#SBATCH -J rsem_alignments
#SBATCH -p standard
#SBATCH -n 1
#SBATCH --cpus-per-task 20
#SBATCH -t 4:00:00
#SBATCH -o ./3-rsem/%a_%A.out
#SBATCH -e ./3-rsem/%a_%A.err

SAMP="opc12-"`printf %02d $SLURM_ARRAY_TASK_ID`

mkdir ./3-rsem/${SAMP}

module load gcc rsem

REF="/home/mds5cg/genomes/grcm38_cDNA_ercc/grcm38_cDNA_rsem"

rsem-calculate-expression -p 20 --bowtie2 --bowtie2-path /home/mds5cg/software/bowtie2/bowtie2-2.3.4.3 --single-cell-prior --paired-end ./2-trim/${SAMP}_R1.trim.fq.gz ./2-trim/${SAMP}_R2.trim.fq.gz $REF ./3-rsem/${SAMP}
```

```bash
for i in {1..96}; do
SAMP="opc12-"`printf %02d $i`
mv ${SAMP}.genes.results ${SAMP}
mv ${SAMP}.isoforms.results ${SAMP}
mv ${SAMP}.transcript.bam ${SAMP}
mv ${SAMP}.stat ${SAMP}
mv ${i}_* ${SAMP}
done
```

# Trimming output
```
SAMPLIST=$(find *log | cut -f 1 -d "." | sort -u)
for SAMP in $SAMPLIST
do
RC=$(grep "Total reads" ${SAMP}.trim.log | cut -f 3 -d " ")
TR=$(grep "Too short" ${SAMP}.trim.log | cut -f 5 -d " ")
RL=$(( $RC - $TR ))
echo -e '|'$SAMP'|'$RC'|'$RL'|'
done

```



All cartridges

Make scratch directory
```bash
cd /scratch/mds5cg
mkdir opc12
cd opc12
mkdir raw
logout
```

Copy files over from first cartridge
```bash
pscp -scp -r "/Volumes/JanesLab/RNA-seq raw data/10162019-Janes_Matt-opc12-142551423/*" mds5cg@rivanna.hpc.virginia.edu:/scratch/mds5cg/opc12/0-raw/



#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH -A janeslab
#SBATCH -p standard
#SBATCH -c 1

SAMPLIST=$(find ../0-raw/*/*/*.fastq.gz | cut -f 5 -d "/" | cut -f 1 -d "_" | sort -u)
for SAMP in $SAMPLIST
do
echo "Started  $SAMP at $(date +%F\ %T)"
R1=$(find ../0-raw/*/*/$SAMP*R1*.fastq.gz | sort)
R2=$(find ../0-raw/*/*/$SAMP*R2*.fastq.gz | sort)
cat $R1 > ./${SAMP}_R1.fastq.gz
cat $R2 > ./${SAMP}_R2.fastq.gz
echo "Finished $SAMP at $(date +%F\ %T)"
done


#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH -A janeslab
#SBATCH -p standard
#SBATCH -c 1
SAMPLIST=$(find ../1-cat/*.fastq.gz | cut -f 3 -d "/" | cut -f 1 -d "_" | sort -u)
for SAMP in $SAMPLIST
do
echo "Started  $SAMP at $(date +%F\ %T)"
~/software/ea-utils/ea-utils-master/clipper/fastq-mcf -q 10 -t 0.01 -k 0 \
    ~/adapters/nextera_adapters.fa \
    ../1-cat/${SAMP}_R1.fastq.gz ../1-cat/${SAMP}_R2.fastq.gz \
    -o ./${SAMP}_R1.trim.fastq.gz -o ./${SAMP}_R2.trim.fastq.gz \
    > ./${SAMP}.trim.log
echo "Finished $SAMP at $(date +%F\ %T)"
done

module load gcc gparallel
find ./1-concat/*fastq.gz | cut -f 3 -d "/" | cut -f 1 -d "_" | sort | parallel -j 12 "~/software/ea-utils/ea-utils-master/clipper/fastq-mcf -q 10 -t 0.01 -k 0 ~/adapters/nextera_adapters.fa ./1-concat/{}_R1.fastq.gz ./1-concat/{}_R2.fastq.gz -o ./2-trim/{}_R1.trim.fq.gz -o ./2-trim/{}_R2.trim.fq.gz > ./2-trim/{}.trim.log"

```


```
#!/bin/bash
# trim_array.sh

#SBATCH -J trim_array
#SBATCH -A janeslab
#SBATCH -p standard
#SBATCH -c 1
#SBATCH -t 1:00:00
#SBATCH -o ./%a_%A.out

SAMP="opc12-"`printf %02d $SLURM_ARRAY_TASK_ID`

echo "Started  $SAMP at $(date +%F\ %T)"
~/software/ea-utils/ea-utils-master/clipper/fastq-mcf -q 10 -t 0.01 -k 0 \
    -o ./${SAMP}_R1.trim.fastq.gz -o ./${SAMP}_R2.trim.fastq.gz
    ~/adapters/nextera_adapters.fa \
    ../1-cat/${SAMP}_R1.fastq.gz ../1-cat/${SAMP}_R2.fastq.gz \
    > ./trim_logs/${SAMP}.trim.log
echo "Finished $SAMP at $(date +%F\ %T)"




REF="/home/mds5cg/genomes/grcm38_cDNA_ercc/grcm38_cDNA_rsem"

rsem-calculate-expression -p 20 --bowtie2 --bowtie2-path /home/mds5cg/software/bowtie2/bowtie2-2.3.4.3 --single-cell-prior --paired-end ./2-trim/${SAMP}_R1.trim.fq.gz ./2-trim/${SAMP}_R2.trim.fq.gz $REF ./3-rsem/${SAMP}





```




```
#!/bin/bash
# rsem_array.sh

#SBATCH -J rsem_array
#SBATCH -p standard
#SBATCH -n 1
#SBATCH --cpus-per-task 20
#SBATCH -t 1:00:00
#SBATCH -o ./slurm_logs/%a_%A.out

SAMP="opc12-"`printf %02d $SLURM_ARRAY_TASK_ID`

mkdir -p ./${SAMP}

module load gcc rsem

REF="/home/mds5cg/genomes/grcm38_cDNA_ercc/grcm38_cDNA_rsem"

rsem-calculate-expression -p 20 --bowtie2 --bowtie2-path /home/mds5cg/software/bowtie2/bowtie2-2.3.4.3 --single-cell-prior --paired-end --output-genome-bam \
../2-trim/${SAMP}_R1.trim.fastq.gz ../2-trim/${SAMP}_R2.trim.fastq.gz \
$REF \
./${SAMP}/${SAMP}



multiqc --interactive
