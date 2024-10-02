SAMPLES = info_df['sample'].tolist()
TRIMMED = info_df.loc[(info_df['experiment'].isin(['ATAC', 'ChIPmentation'])) & (info_df['sequencing_type'] == 'paired-end'), 'sample'].tolist()

seq_dict = [
    {'sample': sample,
    'seq_type': info_df.loc[sample, 'sequencing_type']
    }
    for sample in info_df["sample"]
    ]


def get_fastq(wildcards):
    if wildcards.sample in TRIMMED:
        fastqs = ['%s/_fastq/%s_%s.trimmed.fastq.gz' %(outdir, wildcards.sample, read) for read in ['R1', 'R2']]
    else:
        ## single-end fastqs should not have reads indicator in file names
        reads = ['_R1', '_R2'] if info_df.loc[wildcards.sample, 'sequencing_type'] == 'paired-end' else ['']
        fastqs = ['%s/_fastq/%s%s.fastq.gz' %(outdir, wildcards.sample, read) for read in reads]
    #if info_df.loc[wildcards.sample, 'experiment'] in ['ATAC', 'ChIPmentation']:
    return([fq for fq in fastqs])


def get_targets(wildcards):
    
    bam1 = expand('%s/{sequencing_type}/bam/{sample}.raw.bam' %outdir, zip, 
    sequencing_type = [wc['seq_type'] for wc in seq_dict], 
    sample = [wc['sample'] for wc in seq_dict])
    
    bam2 = expand('%s/{sequencing_type}/bam/{sample}.rmdup.bam' %outdir, zip,
    sequencing_type = [wc['seq_type'] for wc in seq_dict], 
    sample = [wc['sample'] for wc in seq_dict])
    
    bw = expand('%s/{sequencing_type}/bigwig/{sample}.cpm.bw' %outdir, zip,
    sequencing_type = [wc['seq_type'] for wc in seq_dict], 
    sample = [wc['sample'] for wc in seq_dict])
    
    ## include QC logs to rule all to run QC rules 
    ### added to main snakefile to optionally run QC
    
    ### qc = expand('%s/{sequencing_type}/QC/multiqc_log.html' %outdir, sequencing_type = ['single-end', 'paired-end'])
    ### targets = bam1 + bam2 + bw + qc

    targets = bam1 + bam2 + bw
    

    return targets


def get_bt2idx(wildcards):
    build = info_df.loc[wildcards.sample, 'build']
    bt2_idx = os.path.join(bt2_index_dir, build)
    return(bt2_idx)
