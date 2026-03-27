process trim_tag_umi_autodetect{
    tag {meta.id}

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.umitagged.fastq.gz'), emit: umitagged_reads
    tuple val(meta), path('*.num_reads'), emit: num_reads

    script:
    """
    ## - Read1: trim i7 adaptor
    zcat ${reads[0]} > _R1.fastq
    bbduk.sh in=_R1.fastq out=_R1.trm7.fastq literal="${params.Read1Tail},GGGGGGGGGGGGGGGGGGGG" ktrim=r k=20 mink=3 edist=0 ordered=t minlength=0 qtrim=r trimq=10 threads=1 -Xmx4g &> trim.R1.log

    # - prep barcodes using the first 10000 reads, with freq at least 10
    head -n 40000 _R1.trm7.fastq | paste - - - - | cut -f 2 \\
        | cut -c ${params.pos1}-${params.pos2} | grep -ve 'CCCC' -ve 'NNNN' -ve '^\$' | sort | uniq -c \\
        | awk '{if (\$1 > 10) print}' | sort -k1,1nr > _i5seq.freq
    sed 's:.* ::' _i5seq.freq > _barcodes
    # Auto-detect FIXSEQ from data
    effective_fixseq=\$(head -1 _barcodes)
    if [[ -z "\$effective_fixseq" ]]; then
        echo "ERROR: FIXSEQ auto-detection failed -- no barcodes passed filtering at positions ${params.pos1}-${params.pos2}" >&2
        exit 1
    fi
    echo "Auto-detected FIXSEQ: \$effective_fixseq"
    agrep -2 "\$effective_fixseq" _barcodes | head -n 10 > _umi_whitelist

    ##=== Read 2: trim adaptors ===
    zcat ${reads[1]} > _R2.fastq
    bbduk.sh in=_R2.fastq out=_R2.trm5.fastq literal="${params.Read2Tail},GGGGGGGGGGGGGGGGGGGG" ktrim=r k=20 mink=3 edist=0 ordered=t minlength=0 qtrim=r trimq=10 threads=1 -Xmx4g &> trim.R2.log

    #=== trim R2 tail r1seqRC
    r1seqRC=\$(echo \$effective_fixseq | rev | tr 'ATCG' 'TAGC' | sed 's: .*::' | sed 's:\t.*::')
    bbduk.sh in=_R2.trm5.fastq out=_R2.trmr1.fastq restrictright=${params.pos2} literal="\$r1seqRC" ktrim=r k=8 mink=7 edist=0 ordered=t minlength=0 qtrim=r trimq=10 threads=1 -Xmx4g

    ##=== Filter read pairs shorter than barcode+UMI pattern ===
    reformat.sh in=_R1.trm7.fastq in2=_R2.trmr1.fastq out=_R1.filt.fastq out2=_R2.filt.fastq minlength=${params.pos2} threads=1 -Xmx1g

    ##=== UMI extraction + barcode filtering (umi-tools) ===
    umi_tools extract \\
        --bc-pattern='NNNNNNNNNNCCCCCCCC' \\
        --stdin=_R1.filt.fastq \\
        --stdout=${meta.id}.r1.umitagged.fastq \\
        --read2-in=_R2.filt.fastq \\
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
