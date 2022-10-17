#!/bin/bash
#SBATCH -N 1
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=j.m.kaczorowska@amsterdamumc.nl
#SBATCH -t 02:30:00

#copy the required data to scratch
cp -r $HOME/variant_call/output "$TMPDIR"

#Run Lofreq variant caller
for file1 in "$TMPDIR"/output/*marked.bam
do
out=${file1%%.bam}.vcf
lofreq call -f "$TMPDIR"/output/ORF1_up2.fasta -o $out $file1
done

#copy back to home
cp "$TMPDIR"/output/*vcf $HOME/variant_call/output_vcf
