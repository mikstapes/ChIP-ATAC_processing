## alignment w bowtie2 & marking duplicates with samblaster
rule align_bt2:
    input: 
        fastq = get_fastq,
    output: os.path.join(outdir, 'bam/{sample}.bam')
    threads: min(workflow.cores, 12)
    log:
        bowtie2 = os.path.join(outdir, 'logs/{sample}.bowtie2.log'),
        markdup = os.path.join(outdir, 'logs/{sample}.markdup.log')
    params:
        alignment = '--mm --very-sensitive --no-discordant -X 1000',
        bt2_index = get_bt2idx
    shell:
        """
        bowtie2 -x {params.bt2_index} {params.alignment} --threads {threads} \
        -1 {input.fastq[0]} -2 {input.fastq[1]} 2> {log.bowtie2} \
        | samblaster 2> {log.markdup} \
        | samtools view -Sb -o {output} 
        """

## removes PCR dups, apply MAPQ 10, remove chrM, etc & sort + indexing

rule filter_alignments:
    input: os.path.join(outdir, 'bam/{sample}.bam')
    output: 
        bam = os.path.join(outdir, 'bam/{sample}.rmdup.bam'),
        bam_idx = os.path.join(outdir, 'bam/{sample}.rmdup.bam.bai')
    threads: min(workflow.cores, 10)
    params:
        filtering = '-Sb -F 4 -q 10 '
    shell:
        """
        samtools view -h {input} \
        | samblaster --removeDups \
        | grep -v -P '\tchrM\t' \
		| samtools view {params.filtering} \
		| samtools sort -m 2G -@ 5 -T {input}.tmp -o {output.bam};

        samtools index {output.bam} {output.bam_idx} 
        """


## make bigwigs, extend PE reads to fragment length
## extend SE reads to 250 bc most fragments are 200-300bp in size

rule makebw:
    input: os.path.join(outdir, 'bam/{sample}.rmdup.bam')
    output: os.path.join(outdir, 'bigwig/{sample}.cpm.bw')
    threads: min(workflow.cores, 10)
    log: os.path.join(outdir, 'logs/{sample}.bamCoverage.log')
    params: 
        extension = lambda wc: '' if info_df.loc[wc.sample, 'sequencing_type'] == 'paired-end' else '250'
    shell:
        """
        bamCoverage \
        --bam {input} --outFileName {output} \
        --binSize 10 --normalizeUsing CPM --skipNAs --extendReads {params.extension} \
        --numberOfProcessors {threads} 2> {log}
        """
    
