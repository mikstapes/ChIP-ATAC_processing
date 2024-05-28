#TRIMMED = info_df.loc[info_df['experiment'].isin(['ATAC', 'ChIPmentation']), 'sample'].tolist()

rule trim_adapters:
    input: 
        R1 = os.path.join(outdir, '_fastq/{sample}_R1_001.fastq.gz'),
        R2 = os.path.join(outdir, '_fastq/{sample}_R2_001.fastq.gz'),
    output: 
        R1 = os.path.join(outdir, '_fastq/{sample}_R1_001.trimmed.fastq.gz'),
        R2 = os.path.join(outdir, '_fastq/{sample}_R2_001.trimmed.fastq.gz'),
    threads: min(workflow.cores, 8)
    log: os.path.join(outdir, 'logs/{sample}.cutadapt.log')
        #expand('{outdir}/logs/{sample}.cutadapt.log', outdir=outdir, sample = TRIMMED)
        #os.path.join(outdir, 'logs', '{sample}.cutadapt.log')
    params:
        adapt1 = 'CTGTCTCTTATACACATCTCCGAGCCCACGAGAC',
        adapt2 = 'CTGTCTCTTATACACATCTGACGCTGCCGACG',
        cutoff = '--minimum-length 25 --quality-cutoff 20'
    shell:
        """
        cutadapt \
        -a {params.adapt1} -A {params.adapt2} -j {threads} {params.cutoff} \
        -o {output.R1} -p {output.R2} {input.R1} {input.R2} &> {log}
        """
