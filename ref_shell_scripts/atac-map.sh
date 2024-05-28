	#!/bin/bash

set -eo pipefail; echo "ATAC-seq processing shell script "BEGIN at "$(date)"


#### Quick and dirty shell script to pre-process fastq reads from PE ChIP

### Setting variables ####

ADAPT_SEQ1="CTGTCTCTTATACACATCTCCGAGCCCACGAGAC"
ADAPT_SEQ2="CTGTCTCTTATACACATCTGACGCTGCCGACGA"

TOOLS_DIR="/project/MDL_Ibrahim/MP_all/tools"
DIR=$(pwd)

index="/project/MDL_Ibrahim/MP_all/bt2_index/mm10/mm10"

##Executables
export PATH=/home/phan/miniconda3/envs/chip_atac/bin/:$PATH
export PATH="$TOOLS_DIR"/deepTools/bin:$PATH #bamCoverage

PICARD_JAR="java -jar $TOOLS_DIR/picard-2.23.4.jar"

#find ./fastq/ -name '*merged*' -printf '%f\n' | sed -e 's/\.fastq\.gz//' > files
FILES=$(cat npc_files)

### Prepare working directory

for f in $FILES; do
	mkdir -p "$f"/{aln,QC,pipe_out,logs}
done

### Main pre-process pipeline, loops over all files inside fastq/

for f in $FILES; do
	fastqc --outdir "$DIR"/"$f"/QC "$DIR"/fastq/"$f"*.fastq.gz

# -------PE mapping with bowtie2
	cutadapt \
	-a "$ADAPT_SEQ1" -A "$ADAPT_SEQ2" \
	--core 16               \
	--minimum-length 25     \
	--quality-cutoff 20     \
	-o "$DIR"/fastq/"$f"_R1.trimmed.fastq.gz -p "$DIR"/fastq/"$f"_R2.trimmed.fastq.gz \
	"$DIR"/fastq/"$f"_R1.fastq.gz "$DIR"/fastq/"$f"_R2.fastq.gz

# -------PE mapping with bowtie2
	bowtie2 \
	--very-sensitive \
	-X 1000 \
	-p 8 \
	-x "$index" \
	-1 "$DIR"/fastq/"$f"_R1.trimmed.fastq.gz \
	-2 "$DIR"/fastq/"$f"_R2.trimmed.fastq.gz \
	| samtools view -u - \
	| samtools sort -@ 8 -o "$DIR"/"$f"/aln/"$f".sorted.bam

# -------Get mapping/pairing, duplicates & inserts stats. Included dups removal
	samtools flagstat "$DIR"/"$f"/aln/"$f".sorted.bam \
	> "$DIR"/"$f"/QC/"$f".sorted.bam.flagstat \
	2> "$DIR"/"$f"/logs/"$f"_flagstat.log

	$PICARD_JAR MarkDuplicates I="$DIR"/"$f"/aln/"$f".sorted.bam \
	O="$DIR"/"$f"/aln/"$f".noDUP.bam \
	M="$DIR"/"$f"/QC/Picard_dups_metrics.txt \
	REMOVE_DUPLICATES=true \
	2>"$DIR"/"$f"/logs/"$f"_markdups.log

	$PICARD_JAR CollectInsertSizeMetrics \
	I="$DIR"/"$f"/aln/"$f".sorted.bam \
	O="$DIR"/"$f"/QC/"$f"_inserts_QC.txt \
	H="$DIR"/"$f"/QC/"$f"_inserts_plot.pdf \
	M=0.5

# ------Remove chrM reads + add mapping filters (MAPQ10) for final bam
	samtools index \
	"$DIR"/"$f"/aln/"$f".noDUP.bam \
	"$DIR"/"$f"/aln/"$f".noDUP.bam.bai

	samtools view -F 4 -q 10 -h \
	"$DIR"/"$f"/aln/"$f".noDUP.bam \
	| grep -v "chrM" \
	| samtools view -Sb \
	-o "$DIR"/"$f"/pipe_out/"$f".filt.bam


# ------Aggregate all QC stats
	multiqc "$DIR"/"$f"/QC \
	-n "$DIR"/"$f"/QC/"$f"_multqc.html \
	-f -v \
	2> "$DIR"/"$f"/logs/"$f"_multiqc.log

# -------Make bigwigs
	samtools index \
	"$DIR"/"$f"/pipe_out/"$f".filt.bam \
	"$DIR"/"$f"/pipe_out/"$f".filt.bam.bai

	bamCoverage \
	--bam "$DIR"/"$f"/pipe_out/"$f".filt.bam \
	--outFileName "$DIR"/"$f"/pipe_out/"$f".bw \
	--binSize 1 \
	--normalizeUsing CPM \
	--skipNAs \
	--extendReads \
	--ignoreForNormalization chrM \
	--numberOfProcessors 8 \
	2>"$DIR"/"$f"/logs/"$f"_makebw.log

done

exitstat=$?

echo exit status was "$exitstat"

echo "ATAC-seq processing shell script "END at "$(date)"

echo "Pipeline by phan@molgen.mpg.de"
