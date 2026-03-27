#!/bin/bash
# Downloads the hg38 reference genome from UCSC and builds BWA index.
# Run this once before using the pipeline.
#
# Prerequisites: wget (or curl), bwa, samtools
#
# Usage: bash scripts/setup_reference.sh [output_directory]
# Default output: data/reference/

set -euo pipefail

OUTDIR="${1:-data/reference}"
mkdir -p "$OUTDIR"

UCSC_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
FASTA="$OUTDIR/hg38.fa"

if [ -f "$FASTA" ] && [ -f "$FASTA.bwt" ]; then
    echo "Reference genome and BWA index already present."
    ls -lh "$OUTDIR"/hg38.fa*
    exit 0
fi

echo "Downloading hg38 FASTA from UCSC..."
wget -O "$OUTDIR/hg38.fa.gz" "$UCSC_URL"
echo "Decompressing..."
gunzip "$OUTDIR/hg38.fa.gz"

echo "Building FASTA index..."
samtools faidx "$FASTA"

echo "Building BWA index (this may take several hours)..."
bwa index "$FASTA"

echo "Reference genome setup complete: $FASTA"
ls -lh "$OUTDIR"/hg38.fa*
