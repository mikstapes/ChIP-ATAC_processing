rule fastqc:
    input: 
        os.path.join(outdir, '_fq1', '{sample}_R1_001.fastq.gz'),
        os.path.join(outdir, '_fq1', '{sample}_R2_001.fastq.gz')
    output:
        os.path.join(outdir, 'QC', '{sample}_R1_001.fastqc.zip'),
        os.path.join(outdir, 'QC', '{sample}_R2_001.fastqc.zip')
    threads: 2
    message: "Running fastqc on {input}"
    params: 
        outdir: os.path.join(outdir, 'QC')
    shell:
        """
        fastqc -o {params.outdir} --noextract {input[0]} {input[1]}
        """