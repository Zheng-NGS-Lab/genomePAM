process trim_tag_umi{
    tag {meta.id}

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.umitagged.fastq'), emit: umitagged_reads
    tuple val(meta), path('*.num_reads'), emit: num_reads

    script:
    """
    ## - Read1: trim i7 adaptor
    zcat ${reads[0]} > _R1.fastq
    ${params.BBMAPDIR}/bbduk.sh in=_R1.fastq out=_R1.trm7.fastq literal="${params.Read1Tail},GGGGGGGGGGGGGGGGGGGG" ktrim=r k=20 mink=3 edist=0 ordered=t minlength=0 qtrim=r trimq=10 threads=1 -Xmx4g &> trim.R1.log
    ${params.BBMAPDIR}/bbduk.sh in=_R1.trm7.fastq out=_R1.nolinker.fastq forcetrimleft=${params.pos2} minlength=0 threads=1 -Xmx4g

    # - prep barcodes using the first 10000 reads, with freq at least 10
    head -n 40000 _R1.trm7.fastq | paste - - - - | cut -f 2 \
        | cut -c ${params.pos1}-${params.pos2} | grep -ve 'CCCC' -ve 'NNNN' -ve '^\$' | sort | uniq -c \
        | awk '{if (\$1 > 10) print}' | sort -k1,1nr > _i5seq.freq
    sed 's:.* ::' _i5seq.freq > _barcodes
    agrep -2 ${params.FIXSEQ} _barcodes | head -n 10 | sed 's/^/FIXSEQ /' > _2_mismatch

    # get mi
    echo | awk -v pos1=${params.pos1} -v pos2=${params.pos2} -v xNs=${params.xNs} 'NR==FNR {a[\$2]=\$1; next} {
        if (FNR%4==2) { mi = substr(\$1, 1, pos1-1);
            bc = substr(\$1, pos1, ${params.posR2});
            if (length(mi)==0){mi=xNs};
            if (length(bc)==0){bc=xNs};
            if (bc in a) {\$1=":umi:"mi" "a[bc]
        } else{\$1=":umi:"mi" "bc}
    } else {\$0=""};
    print \$0}' _2_mismatch _R1.trm7.fastq > _I2.dx
    sed 1d _I2.dx > _I2

    # - Read1: output i5 barcode Read1
    field_num=`head -1 _R1.nolinker.fastq | awk '{print NF}'`
    if [[ \$field_num -eq 2 ]]; then
        # For Illumina platform
        paste _R1.nolinker.fastq _I2 | tr ' ' '\t' | cut -f 1,3,4 > _R1_I2
    elif [[ \$field_num -eq 1 ]]; then
        # For BGI platform
        paste _R1.nolinker.fastq _I2 | tr -s ' ' | tr ' ' '	' | cut -f 1,2,3 > _R1_I2
    fi
    awk '{if (NF == 3) {print \$0} else {print \$1}}' _R1_I2 > _R1.I2.fastq
    ## trimed Read1
    paste - - - - < _R1.I2.fastq | grep -P "\tFIXSEQ" | cut -f 1,2,4,5,6 | sed 's/\t/:/' | tr '\t' '\n' > _pre.R1.fq


    ##=== Read 2 ===
    zcat ${reads[1]} > _R2.fastq
    ${params.BBMAPDIR}/bbduk.sh in=_R2.fastq out=_R2.trm5.fastq literal="${params.Read2Tail},GGGGGGGGGGGGGGGGGGGG" ktrim=r k=20 mink=3 edist=0 ordered=t minlength=0 qtrim=r trimq=10 threads=1 -Xmx4g &> trim.R2.log

    #=== trim R2 tail r1seqRC
    r1seqRC=\$(echo ${params.FIXSEQ} | rev | tr 'ATCG' 'TAGC' | sed 's: .*::' | sed 's:\t.*::')
    ${params.BBMAPDIR}/bbduk.sh in=_R2.trm5.fastq out=_R2.trmr1.fastq restrictright=${params.pos2} literal="\$r1seqRC" ktrim=r k=8 mink=7 edist=0 ordered=t minlength=0 qtrim=r trimq=10 threads=1 -Xmx4g

    # - Read2: output i5 barcode Read2
    if [[ \$field_num -eq 2 ]]; then
        # For Illumina platform
        paste _R2.trmr1.fastq _I2 | tr ' ' '\t' | cut -f 1,3,4 > _R2_I2
    elif [[ \$field_num -eq 1 ]]; then
        paste _R2.trmr1.fastq _I2 | tr -s ' ' | tr ' ' '	' | cut -f 1,2,3 > _R2_I2
    fi
    awk '{if (NF == 3) {print \$0} else {print \$1}}' _R2_I2 > _R2.I2.fastq
    paste - - - - < _R2.I2.fastq | grep -P "\tFIXSEQ" | cut -f 1,2,4,5,6 | sed 's/\t/:/' | tr '\t' '\n' > _pre.R2.fq


    ## - umitag Read1 and Read 2
    sed "s/::umi:/ 1:N:0:0 /" _pre.R1.fq > _pre.r1.fq.0
    sed "s/::umi:/ 2:N:0:0 /" _pre.R2.fq > _pre.r2.fq.0
    awk '{if(FNR%4==2){\$0=substr(\$1,1,8)}else{\$0=""};print \$0}' _pre.r1.fq.0 | sed 1d > _head1
    awk '{if(FNR%4==2){\$0=substr(\$1,1,8)}else{\$0=""};print \$0}' _pre.r2.fq.0 | sed 1d > _head2
    paste _pre.r1.fq.0 _head1 _head2 | awk '{if(FNR%4==1){\$0=\$1" "\$2" "\$3"_"\$4"_"\$5}else{\$0=\$1};print}' > ${meta.id}.r1.umitagged.fastq
    paste _pre.r2.fq.0 _head1 _head2 | awk '{if(FNR%4==1){\$0=\$1" "\$2" "\$3"_"\$4"_"\$5}else{\$0=\$1};print}' > ${meta.id}.r2.umitagged.fastq


    ## - Remove intermediates to free up space
    rm _*


    ## Check ID consistency
    nR1=\$(wc -l ${meta.id}.r1.umitagged.fastq | sed 's: .*::')
    nR2=\$(wc -l ${meta.id}.r2.umitagged.fastq | sed 's: .*::')
    hidR1=\$(head -n1 ${meta.id}.r1.umitagged.fastq | sed 's: .*::' | sed 's/\\/[12]//')
    hidR2=\$(head -n1 ${meta.id}.r2.umitagged.fastq | sed 's: .*::' | sed 's/\\/[12]//')
    tidR1=\$(tail -n4 ${meta.id}.r1.umitagged.fastq | head -n1 | sed 's: .*::' | sed 's/\\/[12]//')
    tidR2=\$(tail -n4 ${meta.id}.r2.umitagged.fastq | head -n1 | sed 's: .*::' | sed 's/\\/[12]//')

    if [[ \$nR1 == \$nR2 && "\$hidR1" == "\$hidR2" && "\$tidR1" == "\$tidR2" ]]; then
        echo \$(wc -l ${meta.id}.r1.umitagged.fastq | cut -d ' ' -f 1) /4 | bc > ${meta.id}.num_reads
    else
        echo 'Read names of forward and reversed reads do not match, exiting...'
        exit 1
    fi
    """
}
