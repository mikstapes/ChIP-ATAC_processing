## alignment w bowtie2 & marking duplicates with samblaster
rule align_bt2_PE:
    input: fastq = get_fastq
    output: os.path.join(outdir, 'paired-end/bam/{sample}.raw.bam')
    threads: min(workflow.cores, 12)
    log:
        bowtie2 = os.path.join(outdir, 'paired-end/QC/{sample}.bowtie2.log'),
        markdup = os.path.join(outdir, 'paired-end/QC/{sample}.markdup.log')
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

rule align_bt2_SE:
    input: get_fastq
    output: os.path.join(outdir, 'single-end/bam/{sample}.raw.bam')
    threads: min(workflow.cores, 8)
    log:
        bowtie2 = os.path.join(outdir, 'single-end/QC/{sample}.bowtie2.log'),
        markdup = os.path.join(outdir, 'single-end/logs/{sample}.markdup.log')
    params:
        alignment = '--mm --very-sensitive',
        bt2_index = get_bt2idx
    shell:
        """
        bowtie2 -x {params.bt2_index} {params.alignment} --threads {threads} \
        -U {input}  2> {log.bowtie2} \
        | samblaster 2> {log.markdup} \
        | samtools view -Sb -o {output} 
        """

## removes PCR dups, unpaired, remove chrM + mapq 30 -> sort + indexing

rule filter_alignments:
    input: os.path.join(outdir, '{sequencing_type}/bam/{sample}.raw.bam')
    output: 
        bam = os.path.join(outdir, '{sequencing_type}/bam/{sample}.rmdup.bam'),
        bam_idx = os.path.join(outdir, '{sequencing_type}/bam/{sample}.rmdup.bam.bai')
    threads: min(workflow.cores, 10)
    params:
        filtering = lambda wc: '-F 4 -q 30 -f 2' if info_df.loc[wc.sample, 'sequencing_type'] == 'paired-end' else '-Sb -F 4 -q 30'
    shell:
        """
        samtools view -h {input} \
        | samblaster --removeDups \
        | grep -v -P '\tchrM\t' \
		| samtools view -Sb {params.filtering} \
		| samtools sort -m 2G -@ 5 -o {output.bam};

        samtools index {output.bam} {output.bam_idx} 
        """


# ## make bigwigs, extend PE reads to fragment length
# ## extend SE reads to 250 bc most fragments are 200-300bp in size

rule makebw:
    input: os.path.join(outdir, '{sequencing_type}/bam/{sample}.rmdup.bam')
    output: os.path.join(outdir, '{sequencing_type}/bigwig/{sample}.cpm.bw')
    threads: min(workflow.cores, 10)
    log: os.path.join(outdir, '{sequencing_type}/logs/{sample}.bamCoverage.log')
    params: 
        extension = lambda wc: '' if info_df.loc[wc.sample, 'sequencing_type'] == 'paired-end' else '250'
    shell:
        """
        bamCoverage \
        --bam {input} --outFileName {output} \
        --binSize 10 --normalizeUsing CPM --skipNAs --extendReads {params.extension} \
        --numberOfProcessors {threads} 2> {log}
        """
    
