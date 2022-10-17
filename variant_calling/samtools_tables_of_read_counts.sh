#!/bin/bash
#SBATCH -N 1
#SBATCH -t 00:40:00
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=j.m.kaczorowska@amsterdamumc.nl

outdir=~/variant_call/output_samtools
cp -r ~/variant_call/output "$TMPDIR"
cd "$TMPDIR"

### making table with reads
for file2 in "$TMPDIR"/output/*marked.bam
do
output=${file2%%}.tab
samtools idxstats $file2 > $output
done

### Putting the table together
awk '
BEGIN   { FS=OFS="\t" }

FNR==1  { lines[0]=lines[0] OFS FILENAME }

FNR==NR { lines[FNR]=$1 }

        { lines[FNR]=lines[FNR] OFS $3 }

END     { for (i=0;i<=FNR;i++)
              print lines[i]
        }
' "$TMPDIR"/output/*tab > merged_table_samtools.txt

### sort the bam.file
for file in "$TMPDIR"/output/*marked.bam
do
out=${file%%}.sorted.bam
samtools sort $file -o $out
done

### stdout - making the table with statistics
for file3 in "$TMPDIR"/output/*.sorted.bam
do
outtab=${file3%%}.stdout
~/Illumina/bowtie/bbmap/pileup.sh in=$file3 out=$outtab
done

### Copy
cp "$TMPDIR"/output/*.txt $outdir
cp "$TMPDIR"/output/*tab $outdir
cp "$TMPDIR"/output/*stdout $outdir
