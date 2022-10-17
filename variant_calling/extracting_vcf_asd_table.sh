#!/bin/bash
#SBATCH -N 1
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=j.m.kaczorowska@amsterdamumc.nl
#SBATCH -t 00:20:00

module load 2020
module load Java/1.8.0_261
module list

#copy the required data to scratch
cp $HOME/variant_call/output_vcf/round2/*.vcf "$TMPDIR"

for file1 in "$TMPDIR"/*vcf
do
output=${file1%%}.table
~/variant_call/gatk-4.2.5.0/gatk VariantsToTable \
     -V $file1 -F CHROM -F REF -F ALT -F POS -F AF -F DP -F TYPE -GF AD \
     -O $output
done

#copy back
cp "$TMPDIR"/*table $HOME/variant_call/tables
