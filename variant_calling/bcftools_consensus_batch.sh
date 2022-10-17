#!/bin/bash
#SBATCH -N 1
#SBATCH -t 00:30:00
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=j.m.kaczorowska@amsterdamumc.nl

outdir=~/variant_call/output_consensus
cp -r ~/variant_call/output_vcf "$TMPDIR"
cd "$TMPDIR"

### indexing files
for file in *vcf
do
bgzip $file
done

for file in *vcf.gz
do
bcftools index $file
done

### run the bcf consensus
for file in *vcf.gz
do
output=${file%%}.fasta
cat ORF1_up2.fasta | bcftools consensus $file > $output
done

### copy the files back to output
cp *fasta $outdir
