def get_logs(outdir, df):

    ## from https://github.com/tobiaszehnder/process_seq_data/blob/main/Snakefile
    logfiles = []
    for idx, row in df.iterrows():
        prefix = '%s/logs/%s' %(outdir, row['sample'])
        logfiles += [prefix + suffix + '.log' for suffix in ('.cutadapt', '.rmdup', '.bamCoverage')]
        # if row['experiment'] == 'chromatin-accessibility':
        #     logfiles += [prefix + '.full.bam.bowtie2.log']
        # else:
        #     logfiles += [prefix + '.full.Log.final.out']
    return logfiles


def get_targets(df):
    SAMPLES = df['sample'].tolist()
    solexa = expand('%s/_fastq/{sample}_{read}_001.fastq.gz' %outdir, sample = SAMPLES, read = ['R1', 'R2'])
    #fastq = expand('%s/_fastq/{sample}_{read}_001.trimmed.fastq.gz' %outdir, sample = SAMPLES, read = ['R1', 'R2'])
    bam = expand('%s/bam/{sample}.rmdup.bam' %outdir, sample = SAMPLES)
    bamidx = expand('%s/bam/{sample}.rmdup.bam.bai' %outdir, sample = SAMPLES)
    bigwig = expand('%s/bigwig/{sample}.rmdup.cpm.bw' %outdir, sample = SAMPLES)
    logs = get_logs(outdir, df)
    targets = solexa + bam + bamidx + bw +logs

    return targets

def get_fastq(wildcards):
    reads = ['R1', 'R2'] if info_df.loc[wildcards.sample, 'sequencing_type'] == 'paired-end' else ['R1']
    if info_df.loc[wildcards.sample, 'experiment'] in ['ATAC', 'ChIPmentation']:
        fastqs = ['%s/_fastq/%s_%s_001.trimmed.fastq.gz' %(outdir, wildcards.sample, read) for read in reads]
    else:
        fastqs = ['%s/_fastq/%s_%s_001.fastq.gz' %(outdir, wildcards.sample, read) for read in reads]
    
    return([fq for fq in fastqs])


def get_bt2idx(wildcards):
    ref_path = config['ref_dir']
    bt2_idx = os.path.join(ref_path, info_df.loc[wildcards.sample, 'build'])
    return(bt2_idx)

# def parse_sample_df(df, get=('solexa_paths', 'samples')):
#     out = []
#     for _,row in df.iterrows():
#         lib = row['library_number'].split('-')[0]
#         flowildcardsell = row['flow_cell']
#         dir = row['solexa']
#         pattern = re.compile(rf'.*{lib}.*{flowildcardsell}.*\.gz$')
#         if get == 'solexa_paths':
#             out += [os.path.join(dir, f) for f in os.listdir(dir) if re.match(pattern, f)]
#         else:
#             out += ['_'.join(row.iloc[1:5]) + '_' + lib]

#     return(out)


