# Get mapping/pairing, duplicates & inserts stats
rule flagstat:
    input: os.path.join(outdir, '{sequencing_type}/bam/{sample}.raw.bam')
    output: os.path.join(outdir, '{sequencing_type}/QC/{sample}.flagstat')
    threads: 1
    log: os.path.join(outdir, '{sequencing_type}/logs/{sample}.flagstat.log')
    shell:
        """
        samtools flagstat {input} > {output} 2> {log}
        """

rule insertsQC:
    input: os.path.join(outdir, '{sequencing_type}/bam/{sample}.raw.bam')
    output: 
        os.path.join(outdir, '{sequencing_type}/QC/{sample}_inserts.txt'),
        os.path.join(outdir, '{sequencing_type}/QC/{sample}_inserts_plot.pdf')
    threads: 1
    log: os.path.join(outdir, '{sequencing_type}/logs/{sample}.inserts.log')
    shell:
        """
        picard CollectInsertSizeMetrics \
        -I {input} -O {output[0]} -H {output[1]}
        """

# Aggregate all QC stats w multiQC

rule multiQC:
### use output from all QC rules as input for multiQC to initiate them
### output of multiQC rule_all
    input: 
        flagstat = expand('%s/{sequencing_type}/QC/{sample}.flagstat' %outdir, zip,
        sequencing_type = [wc['seq_type'] for wc in seq_dict], 
        sample = [wc['sample'] for wc in seq_dict]),
        inserts = expand('%s/{sequencing_type}/QC/{sample}_inserts.txt' %outdir, zip,
        sequencing_type = [wc['seq_type'] for wc in seq_dict], 
        sample = [wc['sample'] for wc in seq_dict])
    output: os.path.join(outdir, '{sequencing_type}/QC/multiqc_log.html')
    params: 
        QCdir = os.path.join(outdir, '{sequencing_type}/QC')
    threads: 1
    shell:
        """
        multiqc {params.QCdir} -o {params.QCdir} -f -v -n multiqc_log 
        """


# Get sequencing stats
# rule fastqc:
#     ### add 001 to avoid reruning rule for trimmed files
#     ### here sample include both R1/R2 or just R1
#     input: os.path.join(outdir, '_fastq/{sample}_001.fastq.gz')
#     output: os.path.join(outdir, '_fastq/fastQC/{sample}_001.fastqc.zip')
#     threads: 2
#     params: 
#         outdir= os.path.join(outdir, '_fastq/fastqc')
#     shell:
#         """
#         fastqc -o {params.outdir} --noextract {input} 
#         """
