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
        fastqs = ['%s/_fastq/%s_%s_001.trimmed.fastq.gz' %(outdir, wildcards.sample, read) for read in ['R1', 'R2']]
    else:
        reads = ['R1', 'R2'] if info_df.loc[wildcards.sample, 'sequencing_type'] == 'paired-end' else ['R1']
        fastqs = ['%s/_fastq/%s_%s_001.fastq.gz' %(outdir, wildcards.sample, read) for read in reads]
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
    qc = expand('%s/{sequencing_type}/QC/multiqc_log.html' %outdir, sequencing_type = ['single-end', 'paired-end'])
    

    #targets = bam1 + bam2 + bamidx + bw + logs

    targets = bam1 + bam2 + bw + qc
    

    return targets

# sample_mates = [{'sample': sample,
#   'mates': ['R1', 'R2'] if info_df.loc[sample, 'sequencing_type'] == 'paired-end' else ['R1']}
#     for sample in info_df["sample"]
#     ]


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

## convert info into dictionary to parameterization of rules
#samples_dict = df.to_dict(orient = "index")
#SAMPLES = list(samples_dict.keys())

#def getQC(df):
#     ## adapted from https://github.com/tobiaszehnder/process_seq_data
#     logfiles = []
#     log_types = ('.flagstat', '.inserts', 'multiqc')
#     for _, row in df.iterrows():
#         prefix = '%s/%s/logs/%s' %(outdir, row['sample'])
#         logfiles += [prefix + suffix + '.log' for suffix in log_types]
#     return logfiles
