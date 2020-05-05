Aligning raw files from TruSeq stranded mRNA Library Prep

**Trimming adapters**- used Alex's master adapters file (~/adapters/all-illumina-adapters-with-pfos-from-alison.fa)

```bash
find ./concat/*fastq.gz | cut -f 1 -d "_" | cut -f 3 -d "/" | sort | parallel -j 6 "~/software/ea-utils/ea-utils-master/clipper/fastq-mcf -q 10 -t 0.01 -k 0 ~/adapters/all-illumina-adapters-with-pfos-from-alison.fa ./concat/{}_R1.fastq.gz ./concat/{}_R2.fastq.gz -o ./trim/{}_R1.trim.fq.gz -o ./trim/{}_R2.trim.fq.gz > ./trim/{}.trim.log"
```

Saved read numbers before and after trimming by enclosing commands in parentheses and then adding >> TrimmedReads.txt after closed parentheses

Quality check- ran FastQC on .gz files- it worked

**Actual alignment with rsem-calculate-expression**

```bash
rsem-calculate-expression -p 28 --bowtie2 --forward-prob 0 --paired-end ../trim/${SAMP}_R1.trim.fq.gz ../trim/${SAMP}_R2.trim.fq.gz $REF ./${SAMP}/$SAMP
```

^--forward-prob 0 option because we used TruSeq stranded kit and the strand that is amplified is the reverse strand (the strand complimentary to mRNA) so the probability that reads align to the forward strand is zero

**CSV generation**- to copy the .isoform.results files of all samples -
First had to copy .isoform.results from all subdirectories into rsem (parent) directory-

```bash
find . -name '*.isoforms.results' -exec cp -t . {} +
```

Things changed in CSV generation- replaced grcm38 with grch38 and took out # Join mouse genes and ERCCs section of code, also within summarize function changed funs to list (R told me funs was self-deprecated)


*Aligned Ray's to grcm38, aligned mine to grch38 (changed in CSV, etc)
