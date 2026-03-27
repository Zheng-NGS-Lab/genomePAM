process trim_tag_umi{
    tag {meta.id}

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.umitagged.fastq.gz'), emit: umitagged_reads
    tuple val(meta), path('*.num_reads'), emit: num_reads

    script:
    """
    ##=== Paired-end adapter + quality trimming (cutadapt, single pass) ===
    r1seqRC=\$(echo "${params.FIXSEQ}" | rev | tr 'ATCG' 'TAGC')

    cutadapt \\
        -a "${params.Read1Tail}" -a "G{20}" \\
        -A "${params.Read2Tail}" -A "G{20}" -A "\$r1seqRC" \\
        -q 10 \\
        --minimum-length ${params.pos2} \\
        -o _R1.trimmed.fastq -p _R2.trimmed.fastq \\
        ${reads[0]} ${reads[1]} > ${meta.id}.cutadapt.log

    ##=== Barcode whitelist from data ===
    head -n 40000 _R1.trimmed.fastq | paste - - - - | cut -f 2 \\
        | cut -c ${params.pos1}-${params.pos2} | grep -ve 'CCCC' -ve 'NNNN' -ve '^\$' | sort | uniq -c \\
        | awk '{if (\$1 > 10) print}' | sort -k1,1nr > _i5seq.freq
    sed 's:.* ::' _i5seq.freq > _barcodes
    agrep -2 "${params.FIXSEQ}" _barcodes | head -n 10 > _umi_whitelist

    ##=== UMI extraction + barcode filtering (umi-tools) ===
    umi_tools extract \\
        --bc-pattern='NNNNNNNNNNCCCCCCCC' \\
        --stdin=_R1.trimmed.fastq \\
        --stdout=${meta.id}.r1.umitagged.fastq \\
        --read2-in=_R2.trimmed.fastq \\
        --read2-out=${meta.id}.r2.umitagged.fastq \\
        --whitelist=_umi_whitelist \\
        --log=${meta.id}.umi_extract.log

    ##=== Count + compress + cleanup ===
    echo \$(wc -l ${meta.id}.r1.umitagged.fastq | cut -d ' ' -f 1) /4 | bc > ${meta.id}.num_reads
    gzip ${meta.id}.r1.umitagged.fastq
    gzip ${meta.id}.r2.umitagged.fastq
    rm _*
    """
}
