## alignment w bowtie2 & marking duplicates with samblaster
#@ WIP: take indexing output as input to prompt indexing before aln
rule align_bt2_PE:
    input: 
        fastq = get_fastq,
        index = get_bt2idx
    output: os.path.join(outdir, 'paired-end/bam/{sample}.raw.bam')
    threads: min(workflow.cores, 8)
    log:
        bowtie2 = os.path.join(outdir, 'paired-end/QC/{sample}.bowtie2.log'),
        markdup = os.path.join(outdir, 'paired-end/QC/{sample}.markdup.log')
    params:
        alignment = config['bt2_params'],
        fragment_size = config['PE_fragment_size'],
        base_idx = lambda wc: info_df.loc[wc.sample, 'build']
    shell:
        """
        bowtie2 -x {input.index}/{params.base_idx} \
        {params.alignment} -X {params.fragment_size} --no-discordant \
        --threads {threads} \
        -1 {input.fastq[0]} -2 {input.fastq[1]} 2> {log.bowtie2} \
        | samblaster 2> {log.markdup} | samtools view -Sb -o {output} 
        """

rule align_bt2_SE:
    input: 
        fastq = get_fastq,
        index = get_bt2idx
    output: os.path.join(outdir, 'single-end/bam/{sample}.raw.bam')
    threads: min(workflow.cores, 4)
    log:
        bowtie2 = os.path.join(outdir, 'single-end/QC/{sample}.bowtie2.log'),
        markdup = os.path.join(outdir, 'single-end/logs/{sample}.markdup.log')
    params:
        alignment = config['bt2_params'],
        base_idx = lambda wc: info_df.loc[wc.sample, 'build']
    shell:
        """
        bowtie2 -x {input.index}/{params.base_idx} \
        {params.alignment} --threads {threads} \
        -U {input.fastq}  2> {log.bowtie2} \
        | samblaster --ignoreUnmated 2> {log.markdup} | samtools view -Sb -o {output} 
        """

## removes PCR dups, unpaired, remove chrM + mapq 30 -> sort + indexing

rule filter_alignments:
    input: os.path.join(outdir, '{sequencing_type}/bam/{sample}.raw.bam')
    output: 
        bam = os.path.join(outdir, '{sequencing_type}/bam/{sample}.rmdup.bam'),
        bam_idx = os.path.join(outdir, '{sequencing_type}/bam/{sample}.rmdup.bam.bai')
    threads: min(workflow.cores, 10)
    params:
        mapq = config['MAPQ'], 
        rmdups = lambda wc: '--removeDups --ignoreUnmated' if info_df.loc[wc.sample, 'sequencing_type'] == 'single-end' else '--removeDups',
        filtering = lambda wc: '-F 4 -f 2' if info_df.loc[wc.sample, 'sequencing_type'] == 'paired-end' else '-Sb -F 4'
    shell:
        """
        samtools view -h {input} \
        | samblaster {params.rmdups} \
        | grep -v -P '\tchrM\t' \
        | samtools view -Sb {params.filtering} -q {params.mapq} \
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
        extension = lambda wc: '' if info_df.loc[wc.sample, 'sequencing_type'] == 'paired-end' else '250',
        bw_bin = config['bw_binsize'],
        bw_norm = config['bw_norm']
    shell:
        """
        bamCoverage \
        --bam {input} --outFileName {output} \
        --binSize {params.bw_bin} \
        --normalizeUsing {params.bw_norm} \
        --extendReads {params.extension} \
        --skipNAs --numberOfProcessors {threads} 2> {log}
        """
    
