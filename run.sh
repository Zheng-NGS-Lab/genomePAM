#!/bin/bash
# Run the genomePAM pipeline locally with Docker.
#
# Prerequisites:
#   1. Docker installed and running
#   2. Nextflow installed
#   3. SRA Toolkit (fasterq-dump) for downloading example data
#
# Usage: bash run.sh [parameters.yml]

set -euo pipefail

PARAMS_FILE="${1:-parameters.yml}"
IMAGE_NAME="genomepam:latest"
PROJECT_DIR="$(cd "$(dirname "$0")" && pwd)"

# Build Docker image if not present
if ! docker image inspect "$IMAGE_NAME" &>/dev/null; then
    echo "Building Docker image..."
    docker build -t "$IMAGE_NAME" .
fi

# Download sample data from NCBI SRA if not present
if [ ! -d "data/fastq" ] || [ -z "$(ls data/fastq/*.fastq.gz 2>/dev/null)" ]; then
    echo "Downloading sample data from NCBI SRA (SRR33421097)..."
    mkdir -p data/fastq
    fasterq-dump --split-3 --outdir data/fastq SRR33421097
    gzip data/fastq/SRR33421097_1.fastq data/fastq/SRR33421097_2.fastq
    mv data/fastq/SRR33421097_1.fastq.gz data/fastq/SRR33421097_R1.fastq.gz
    mv data/fastq/SRR33421097_2.fastq.gz data/fastq/SRR33421097_R2.fastq.gz
fi

# Download reference genome if not present
if [ ! -f "data/reference/hg38.fa" ] || [ ! -f "data/reference/hg38.fa.bwt" ]; then
    echo "Downloading reference genome..."
    bash scripts/setup_reference.sh
fi

# Generate parameters with absolute paths for Docker compatibility
echo "Generating parameters with absolute paths..."
sed \
    -e "s|FQDIR:.*|FQDIR: \"${PROJECT_DIR}/data/fastq\"|" \
    -e "s|OUTDIR:.*|OUTDIR: \"${PROJECT_DIR}/results\"|" \
    -e "s|genome_fasta:.*|genome_fasta: \"${PROJECT_DIR}/data/reference/hg38.fa\"|" \
    "$PARAMS_FILE" > parameters_local.yml

# Create output directory
mkdir -p results

echo "Starting genomePAM pipeline..."
nextflow run main.nf \
    -params-file parameters_local.yml \
    -profile local \
    -with-report results/run_report.html \
    -with-trace results/run_trace.txt \
    -resume
