# ChIP-ATAC_processing

Snakemake pipeline for processing ChIP-seq and ATAC-seq data

---

## Description

Processing pipeline  for ChIP and/or ATAC-seq data written as a Snakemake workflow. 
- Adapters are trimmed from fastq reads before single-end or paired-end alignment with bowtie2.
- Aligned bam files are filtered to remove duplicates, chrM reads and other low-quality reads before generating bigwigs. 
- Quality-control checks at every step of the workflow

## Quick start

The pipeline can be used when provided the following inputs:
1. A `config.yaml` file containing all global variables/paths. See example in `config/config.yaml`
2. A `.tsv` file containing all sample metadata using the template provided [here](https://docs.google.com/spreadsheets/d/e/2PACX-1vQxLXG3ihSa-RvRP5OD1GscYrD-1psD76lg0kdIpl_LsCbcINwwTvYKCee3b7xBlkrnshiWcnYywvk6/pub?gid=0&single=true&output=tsv) 

For example, running the pipeline with default parameters using 10 cores:
` ./process_chipatac -c config/config.yaml -p 10 `


## Usage
Run `./process_chipatac -h ` for detailed usage:

```
Usage: process_chipatac -c CONFIG_FILE (optional: -p THREADS -n DRY_RUN)

Processing pipeline for ChIP & ATAC-seq data from fastq to bigwig.
Author: Mikie Phan (augenlidlos@proton.me)

OPTIONS:
   -c CONFIG    Path to config file in yaml format.
   -p THREADS   Number of threads. Default = 1.
   -s SOLEXA    Flag to indicate data generated at MPI-MG SeqCore.
                For other data sources (default), store fastqs in _fastq/ within output directory (as indicated in <config>.yaml).
   -q QC        Run QC for ChIP/ATAC processing (0=no QC, 1=run QC). Default = 1.
   -n DRY_RUN   Flag to do a dry-run, good to test the pipeline before executing jobs
   -! SMK_ARGS  Optional arguments for Snakemake in quotes, i.e. -! "--debug-dag --unlock -p"
```

- MPI-MG users will need to add the `-s` flag to indicate Solexa file paths
- Changes to mapping/filtering parameters can be done within the provided `config.yaml`

