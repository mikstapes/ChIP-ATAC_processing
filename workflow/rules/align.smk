rule bt2_paired_end:
    input: get_fastq(info_df, aligner = 'bowtie2')
    output: 
        bam = expand('%s/bam/{sample}.sorted.bam' %outdir, sample = SAMPLES),
        #bam_idx = expand('%s/bam/{sample}.rmdup.bam.bai' %outdir, sample = SAMPLES)
    threads: min(workflow.cores, 12)
    params:
        bt2_params = '--mm --very-sensitive --no-discordant -X 1000',
        bt2_index = get_bt2idx(info_df)
    shell:
        """
        bowtie2 -x {params.bt2_index} {params.bt2_params} --threads {threads} \
        -1 {input[0]} -2 {input[1]} \
        | samtools -u - | samtools sort -o {output.bam}
        """
