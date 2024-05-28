#!/bin/bash

set -eo pipefail; echo "ChIPmentation processing shell script "BEGIN at "$(date)"

### Setting variables ####

ADAPT_SEQ1="CTGTCTCTTATACACATCTCCGAGCCCACGAGAC"
ADAPT_SEQ2="CTGTCTCTTATACACATCTGACGCTGCCGACGA"

TOOLS_DIR="/project/MDL_Ibrahim/MP_all/tools"
DIR=$(pwd)

index="/project/MDL_Ibrahim/ATAC/index/gg6/galGal6"

##Executables & variables def
export PATH=/home/phan/miniconda3/envs/chip_atac/bin/:$PATH

export PATH="$TOOLS_DIR"/deepTools/bin:$PATH #bamCoverage

PICARD_JAR="java -jar $TOOLS_DIR/picard-2.23.4.jar"

FILES=$(cat files)


### Prepare working directory

if [ ! -d ]


for f in $FILES; do
	mkdir -p "$f"/{aln,QC,logs}
done

### Main pre-process pipeline, loops over all files inside fastq/

for f in $FILES; do
	fastqc --outdir "$f"/QC fastq/"$f"*.fastq.gz

	cutadapt \
	-a "$ADAPT_SEQ1" -A "$ADAPT_SEQ2" \
	--core 16               \
	--minimum-length 20     \
	--quality-cutoff 20     \
	-o fastq/"$f"_R1.trimmed.fastq.gz -p fastq/"$f"_R2.trimmed.fastq.gz \
	fastq/"$f"_R1_001.fastq.gz fastq/"$f"_R2_001.fastq.gz

# -------PE mapping with bowtie2
	bowtie2 \
	--very-sensitive \
	-X 700 \
	-p 8 \
	-x "$index" \
	-1 fastq/"$f"_R1.trimmed.fastq.gz \
	-2 fastq/"$f"_R2.trimmed.fastq.gz \
	| samtools view -u - \
	| samtools sort -@ 8 - > "$f"/aln/"$f".sorted.bam

# -------Get mapping/pairing & duplicate stats. Dups removal

	samtools flagstat "$f"/aln/"$f".sorted.bam \
	> "$f"/QC/"$f".sorted.bam.flagstat \
	2> "$f"/logs/"$f"_flagstat.log

	$PICARD_JAR MarkDuplicates I="$f"/aln/"$f".sorted.bam \
	O="$f"/aln/"$f".noDUP.bam \
	M="$f"/Picard_reports/Picard_dups_metrics.txt \
	REMOVE_DUPLICATES=true \
	2>"$f"/logs/"$f"_markdups.log

# ------Remove chrM reads + add mapping filters (MAPQ10) for final bam
	samtools index \
	"$f"/aln/"$f".noDUP.bam \
	"$f"/aln/"$f".noDUP.bam.bai

	samtools view -F 4 -q 10 -h \
	"$f"/aln/"$f".noDUP.bam \
	| grep -v "chrM" \
	| samtools view -Sb \
	-o "$f"/pipe_out/"$f".filt.bam

# ------Aggregate all QC stats
	multiqc "$f"/QC \
	-o "$f"/QC/"$f"_multqc.html \
	-f -v \
	2> "$f"/logs/"$f"_multiqc.log

# -------Make bigwigs
	samtools index \
	"$f"/pipe_out/"$f".filt.bam \
	"$f"/pipe_out/"$f".filt.bam.bai

	bamCoverage \
	  --bam "$f"/pipe_out/"$f".filt.bam \
	  --outFileName "$f"/pipe_out/"$f".filt.bw \
	  --binSize 10 \
	  --normalizeUsing RPKM \
	  --skipNAs \
	  --extendReads \
	  --ignoreForNormalization chrM \
	  --numberOfProcessors 8 \
	  2>"$DIR"/"$f"/logs/"$f"_makebw.log

done


exitstat=$?

echo "exit status was $exitstat"

echo "ChIPmentation processing shell script "END at "$(date)"

echo "Pipeline by phan@molgen.mpg.de"
