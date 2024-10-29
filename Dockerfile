# hash:sha256:455cbd41a0b28348f392b6a7a868b0cdfa5fe06bfcc79e4cbe9b3abc0c2eb501
FROM registry.codeocean.com/codeocean/ubuntu:20.04.2

RUN apt-get update && apt-get install -y \
    wget \
    git \
    unzip \
    build-essential \
    bzip2

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
    bash miniconda.sh -b -p $HOME/miniconda && \
    rm miniconda.sh
# Add Miniconda to the PATH
ENV PATH="/root/miniconda/bin:${PATH}"

# Update Conda and install Nextflow
RUN conda update -n base -c defaults conda

# Install Nextflow
RUN conda install -c bioconda nextflow

# Install GUIDE-seq
RUN git clone https://github.com/tsailabSJ/guideseq.git

# Install BBMap
RUN wget https://sourceforge.net/projects/bbmap/files/BBMap_39.10.tar.gz && \
    tar -xzvf BBMap_39.10.tar.gz && \
    rm BBMap_39.10.tar.gz && \
    mv bbmap /usr/local/bin

# Install snpEff
RUN wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip && \
    unzip snpEff_latest_core.zip && \
    rm snpEff_latest_core.zip

# Install BWA
RUN conda install bioconda::bwa

# Set the default command to run when the container starts
CMD ["/bin/bash"]
ARG DEBIAN_FRONTEND=noninteractive
