# genomePAM
Repository for genomePAM

## Setup
Please install the following programs
## Pre-requisites
- Nextflow (https://www.nextflow.io/)
- GUIDE-seq (https://github.com/tsailabSJ/guideseq)
- BBMap (https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmap-guide/)
- snpEff (http://pcingola.github.io/SnpEff/)
- BWA (https://github.com/lh3/bwa)
- BWA-indexed reference genome (tested on GRCh37)
- conda


## Setup
Please install the required programs and conda environment

Conda environment YAML can be found in 
```conda/genomePAM.yml```

Please use BWA to index the reference genome
```bwa index hg38.fna```

### Config
Please change the path to BBMap, guideseq, snpEff and bwa in `nextflow.config`
```
guide_seq_conda = "path/to/guideseq/env"
BBMAPDIR = "/path/to/bbmap"
GUIDESEQDIR = "/path/to/guideseq"
SNPEFFDIR = "/path/to/snpEff"
BWA = "/path/to/bwa"
```
Set the bwa indexed hg38's path
```
// Reference genomes (bwa indexed)
hg38 = "/path/to/bwa_indexed/hg38.fna"
```

## Usage
### Inputs
Path to input directories and corresponding parameters has to be specified in a `parameters.yml` file:

1. `FQDIR`: Path to input directories containing demultiplexed pair-end FASTQ
2. `OUTDIR`: Path to output directories
3. `BWATHREADS`: Number of threads used in the BWA alignment step
4. `Read1Tail`: Custom sequence added to the tail of read 1
5. `Read2Tail`: Custom sequence added to the tail of read 2
6. `pos1`: 
7. `pos2`: 
8. `posR2`: 
9. `xNs`: 
10. `FIXSEQ`: 
11. `GENOME`: Path to a BWA-indexed reference genome
12. `HGVER`: hg19 or hg38 (only tested on hg19)

### Outputs
1. BWA alignment in BAM
2. Table of identified offtarget sites (raw and annotated)
3. Visualization of identified offtargets and PAM sequence logo
4. MultiQC reports of raw FASTQ and trimmed+consolidated FASTQ

### Command
```
nextflow run main.nf -params-file parameters.yml -with-report run_report.html
```