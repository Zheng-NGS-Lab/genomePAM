process align_identify{
    tag {meta.id}

    conda params.guide_seq_conda

    input:
    tuple val(meta), path(consolidated_reads)

    output:
    tuple val(meta), path('*.bam'), emit: aligned_bam
    tuple val(meta), path('*.bam.bai'), emit: aligned_bam_idx
    tuple val(meta), path('*_identifiedOfftargets.txt'), emit: identified_offtargets

    beforeScript "mkdir aligned"
    afterScript "rmdir identified aligned"

    script:
    """
    # Add GUIDEseq env bin to PATH
    PATH=${params.guide_seq_conda}/bin:\$PATH
    
    #####Main
    TargetSeq=\$(echo "${params.AssaySpec}"| tr '[:lower:]' '[:upper:]'|sed 's/_//g')
    
    reference=${params.hg38}

    # Align the consolidated FASTQ to reference
    ${params.BWA} mem -t ${params.BWATHREADS} \$reference \\
        ${consolidated_reads[0]} ${consolidated_reads[1]} > aligned/${meta.id}.sam
    
    # Identify off-target sites with guideseq.py
    python3 ${params.GUIDESEQDIR}/guideseq.py identify \\
        --aligned ${meta.id}.sam \\
        --target_sequence \$TargetSeq --genome \$reference --outfolder .
    
    # Extract results to current directory
    mv ./identified/* ./aligned/* .

    # sam to sorted bam and index
    samtools view -t ${params.BWATHREADS} -bS ${meta.id}.sam > _.bam
    samtools sort --threads ${params.BWATHREADS} _.bam -o ${meta.id}.bam
    samtools index ${meta.id}.bam
    
    # Clean up
    rm -rf ${meta.id}.sam _.bam
    """
}
