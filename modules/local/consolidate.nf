process consolidate{
    tag {meta.id}

    input:
    tuple val(meta), path(umitagged_reads)

    output:
    tuple val(meta), path('*.consolidated.fastq.gz'), emit: consolidated_reads

    script:
    """
    paste - - - - < ${umitagged_reads[0]} | tr ' ' '\t' | gawk '{print \$3,\$6,\$2,\$4,\$5}' > __r1
    LC_ALL=C sort __r1 > __r1srt
    gawk '{if(\$1==pre1){uniq++}else{print uniq,pre0;uniq=1}; pre1=\$1; pre0=\$0;}' __r1srt > __r1srtCnt
    sed 1d __r1srtCnt | awk '{OFS="\t";print "@"\$2"_"\$1" "\$4,\$5,\$6,\$3}' | tr '\t' '\n' > ${meta.id}.r1.consolidated.fastq

    paste - - - - < ${umitagged_reads[1]} | tr ' ' '\t' | gawk '{print \$3,\$6,\$2,\$4,\$5}' > __r2
    LC_ALL=C sort __r2 > __r2srt
    gawk '{if(\$1==pre1){uniq++}else{print uniq,pre0;uniq=1}; pre1=\$1; pre0=\$0;}' __r2srt > __r2srtCnt
    sed 1d __r2srtCnt | awk '{OFS="\t";print "@"\$2"_"\$1" "\$4,\$5,\$6,\$3}' | tr '\t' '\n' > ${meta.id}.r2.consolidated.fastq

    # gzip the consolidated FASTQ
    gzip *.consolidated.fastq

    # Clean-up
    rm __r1* __r2*
    """
}