#!/bin/bash
#SBATCH -N 1
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=j.m.kaczorowska@amsterdamumc.nl
#SBATCH -t 06:00:00

### Load modules; GATK requires downgraded version of Java.

module load 2020
module load Java/1.8.0_261
module list

### Copy the required data to scratch
## input are the read files, trimmed with trimmomatic in fastq format.
cp -r $HOME/variant_call/input2 "$TMPDIR"

### 1 Alignment and sam preprocessing
### 1.1 Run BWA alignment
### for now I use the standard settings of BWA; change if necessary.

for file1 in "$TMPDIR"/input2/*1P.fastq.gz
do
file2=${file1/1P/2P}
out=${file1%%.fastq.gz}_aligned

bwa mem "$TMPDIR"/input2/ORF1 $file1 $file2 > $out
done

### 1.2 Perform sam to bam conversion
for i in "$TMPDIR"/input2/*_aligned
do
out=${i%%.sam}conv.bam

samtools view -bS $i > $out
done

### 1.3 Add or replace the read groups

for j in "$TMPDIR"/input2/*conv.bam
do
out=${j%%conv.bam}_rg.bam

$HOME/variant_call/gatk-4.2.5.0/gatk AddOrReplaceReadGroups \
       I=$j\
       O=$out \
       RGID=4 \
       RGLB=lib1 \
       RGPL=ILLUMINA \
       RGPU=unit1 \
       RGSM=20
done

### 2.1 - mark duplicates in the sam file

for k in "$TMPDIR"/input2/*rg.bam
do
out=${k%%_rg.bam}_marked.bam
$HOME/variant_call/gatk-4.2.5.0/gatk MarkDuplicatesSpark \
 -I $k \
 -O $out
done

### LAST PART copy the files back to home - only marked bam and bai files

#give me back my files
cp "$TMPDIR"/input2/*marked* $HOME/variant_call/output
