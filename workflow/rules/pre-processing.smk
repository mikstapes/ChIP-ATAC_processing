rule bowtie2_index:
    input: os.path.join(config['references'], 'fasta', '{build}.fa')
    output: directory('%s/{build}' %bt2_index_dir)
    shell:
        """
        [[ -d {output} ]] || mkdir -p {output}
        bowtie2-build {input} {output}/{wildcards.build}
        """

for fname in TRIMMED:
    ## no name so can re-run the same rule
    rule:
        input: 
            R1 = os.path.join(outdir, '_fastq/%s_R1.fastq.gz' %fname),
            R2 = os.path.join(outdir, '_fastq/%s_R2.fastq.gz' %fname),
        output: 
            R1 = os.path.join(outdir, '_fastq/%s_R1.trimmed.fastq.gz' %fname),
            R2 = os.path.join(outdir, '_fastq/%s_R2.trimmed.fastq.gz' %fname),
        threads: min(workflow.cores, 8)
        log: os.path.join(outdir, 'logs/%s.cutadapt.log' %fname)
            #expand('{outdir}/logs/%s.cutadapt.log', outdir=outdir, sample = TRIMMED)
            #os.path.join(outdir, 'logs', '%s.cutadapt.log')
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