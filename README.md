![genomePAM](resources/img/genomePAM_logo.png)

# GenomePAM

Nextflow DSL2 pipeline for genome-wide PAM discovery via GUIDE-seq off-target analysis.

All dependencies are containerized in a single Docker image -- no Conda environments, manual tool installations, or path configuration required.

## Prerequisites

- [Docker](https://www.docker.com/) (running)
- [Nextflow](https://www.nextflow.io/) (>= 21.10)
- [SRA Toolkit](https://github.com/ncbi/sra-tools) (for downloading example data from NCBI SRA)

## Quick Start

```bash
# Clone the repo
git clone https://github.com/Zheng-NGS-Lab/genomePAM.git
cd genomePAM

# Run the pipeline (builds Docker image, downloads data, runs everything)
bash run.sh
```

`run.sh` handles the full workflow:
1. Builds the Docker image (`genomepam:latest`) if not already present
2. Downloads sample FASTQ from NCBI SRA (accession [SRR33421097](https://www.ncbi.nlm.nih.gov/sra/?term=SRR33421097))
3. Downloads and BWA-indexes the hg38 reference genome
4. Generates `parameters_local.yml` with absolute paths
5. Runs the Nextflow pipeline with `-resume` support

## Manual Execution

If you prefer to run steps individually:

### 1. Build the Docker image

```bash
docker build -t genomepam:latest .
```

### 2. Prepare data

Place paired-end FASTQ files in `data/fastq/` with the naming convention `*_R1*.fastq.gz` and `*_R2*.fastq.gz`.

Place a BWA-indexed hg38 reference in `data/reference/`:
```bash
bash scripts/setup_reference.sh
# Or manually:
# bwa index hg38.fa
# samtools faidx hg38.fa
```

### 3. Configure parameters

Edit `parameters.yml` to match your experiment:

| Parameter | Description |
|-----------|-------------|
| `FQDIR` | Path to directory containing paired-end FASTQ files |
| `OUTDIR` | Path to output directory |
| `GENOME` | Genome identifier (`hg38` enables annotation via snpEff) |
| `genome_fasta` | Path to BWA-indexed reference FASTA |
| `AssaySpec` | Target sequence and PAM: `SPACER_NNNN` (3' PAM) or `NNNN_SPACER` (5' PAM) |
| `BWATHREADS` | Number of BWA alignment threads (default: 4) |
| `PAMlen` | PAM length for GenomePAM analysis (default: 4) |
| `PAMpos` | PAM position: 3 for 3' PAM, 5 for 5' PAM (default: 3) |

Trimming parameters (`Read1Tail`, `Read2Tail`, `pos1`, `pos2`, `posR2`, `xNs`) have sensible defaults and typically do not need modification. `FIXSEQ` defaults to `"auto"`, which detects the conserved barcode sequence from the data at R1 positions `pos1`-`pos2`. Set `FIXSEQ` to an explicit 8-mer in `parameters.yml` to override auto-detection.

### 4. Run the pipeline

```bash
nextflow run main.nf \
    -params-file parameters_local.yml \
    -profile local \
    -with-report results/run_report.html \
    -with-trace results/run_trace.txt \
    -resume
```

## Pipeline Processes (11 total)

| Process | Description |
|---------|-------------|
| FASTQC_PRE | QC on raw sequencing reads |
| MULTIQC_PRE | Aggregate pre-trim QC report |
| trim_tag_umi | Adapter trimming (cutadapt), UMI extraction + barcode filtering (umi-tools) |
| trim_tag_umi_autodetect | Same as trim_tag_umi but auto-detects FIXSEQ barcode from data (BBDuk) |
| FASTQC_POST | QC on UMI-extracted reads |
| MULTIQC_POST | Aggregate post-trim QC report |
| align_identify | BWA alignment, umi-tools dedup, off-target identification |
| annotate | Annotate off-target sites with snpEff (hg38 only) |
| svg_visualize | SVG visualization of off-target sites |
| genomePAM | PAM frequency analysis and reporting |
| visualize | PAM sequence logos, cumulative read count plots, 4-position heatmaps |

## Outputs

Results are written to the `OUTDIR` directory (default: `results/`):

```
results/
  align_identify/     # BAM alignment, identified off-target sites
  annotate/           # snpEff-annotated off-target table
  svg_visualize/      # SVG off-target visualization
  genomePAM/          # PAM discovery tables and HTML report
  visualize/          # PAM sequence logos, heatmaps, stats CSV
  qc/                 # Pre- and post-trim MultiQC reports
  run_report.html     # Nextflow execution report
  run_trace.txt       # Nextflow process trace
```

Key outputs for PAM analysis:
- `genomePAM/*_GenomePAM_Tab.html` -- GenomePAM report (PAM frequency table)
- `visualize/*_visualize.pdf` -- Sequence logos and cumulative read count plots
- `visualize/*_PAM_PM_1-4.pdf` -- 4-position PAM heatmap (perfect matches)
- `visualize/*_PAM_MM_1-4.pdf` -- 4-position PAM heatmap (mismatches)
- `visualize/allLib_stats.csv` -- Per-library summary statistics

## Changelog

See [CHANGELOG.md](CHANGELOG.md) for a detailed history of changes.

## License

See [LICENSE](LICENSE).
