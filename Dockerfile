FROM condaforge/mambaforge:latest

LABEL maintainer="genomePAM pipeline"
LABEL description="All-in-one container for the genomePAM CRISPR PAM discovery pipeline"

# System dependencies: agrep (approximate grep), gawk, bc, procps (for ps command)
RUN apt-get update && apt-get install -y --no-install-recommends \
    tre-agrep \
    gawk \
    bc \
    procps \
    && rm -rf /var/lib/apt/lists/* \
    && ln -sf /usr/bin/tre-agrep /usr/bin/agrep

# Bioinformatics tools, Python, and R in a single solve
RUN mamba install -y -c bioconda -c conda-forge -c defaults \
    bwa \
    samtools \
    bbmap=39.10 \
    fastqc=0.12.1 \
    snpeff \
    cutadapt \
    umi_tools \
    python=3.10 \
    biopython \
    pyfaidx \
    pysam \
    htseq \
    regex \
    svgwrite \
    pandas \
    numpy \
    scipy \
    pyyaml \
    r-base=4.3 \
    r-ggplot2 \
    r-ggseqlogo \
    r-patchwork \
    r-plyr \
    r-gt \
    r-dplyr \
    r-stringr \
    r-readr \
    r-purrr \
    r-glue \
    pip \
    && pip install --no-cache-dir multiqc \
    && mamba clean -a -y

# Download snpEff hg38 database (needed for annotation step)
RUN snpEff download -v hg38

WORKDIR /app
CMD ["/bin/bash"]
