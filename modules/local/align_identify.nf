process align_identify{
    tag {meta.id}

    input:
    tuple val(meta), path(umitagged_reads)
    path(genome_fasta)
    path(genome_index)

    output:
    tuple val(meta), path('*.bam'), emit: aligned_bam
    tuple val(meta), path('*.bam.bai'), emit: aligned_bam_idx
    tuple val(meta), path('*_identifiedOfftargets.txt'), emit: identified_offtargets

    beforeScript "mkdir -p aligned identified"
    afterScript "rm -rf identified aligned"

    script:
    """
    TargetSeq=\$(echo "${params.AssaySpec}" | tr '[:lower:]' '[:upper:]' | sed 's/_//g')

    # Align the UMI-tagged FASTQ to reference (index files staged alongside fasta)
    bwa mem -t ${task.cpus} ${genome_fasta} \
        ${umitagged_reads[0]} ${umitagged_reads[1]} > aligned/${meta.id}_aligned.sam

    # Convert to BAM and sort by coordinate
    samtools view -@ ${task.cpus} -bS aligned/${meta.id}_aligned.sam | \
        samtools sort --threads ${task.cpus} -o aligned/${meta.id}_sorted.bam

    # Index sorted BAM for umi-tools
    samtools index aligned/${meta.id}_sorted.bam

    # Deduplicate with umi-tools using paired-end mode
    umi_tools dedup \
        --paired \
        --stdin=aligned/${meta.id}_sorted.bam \
        --stdout=aligned/${meta.id}_dedup.bam \
        --output-stats=${meta.id}_umi_stats

    # Convert deduplicated BAM back to SAM for guideseq.py
    samtools view -@ ${task.cpus} -h aligned/${meta.id}_dedup.bam > aligned/${meta.id}.sam

    # Identify off-target sites with guideseq.py
    python3 ${params.GUIDESEQDIR}/guideseq.py identify \
        --aligned aligned/${meta.id}.sam \
        --target_sequence \$TargetSeq --genome ${genome_fasta} --outfolder .

    # Extract results to current directory
    mv ./identified/* .

    # Move deduplicated BAM to output and index
    mv aligned/${meta.id}_dedup.bam ${meta.id}.bam
    samtools index ${meta.id}.bam
    """
}
