process genomePAM{
    tag{meta.id}

    errorStrategy 'ignore'

    input:
    tuple val(meta), path(identified_offtargets)

    output:
    tuple val(meta), path('*_GenomePAM_Tab.html'), emit: genomePAM_html
    tuple val(meta), path('*_GenomePAM.txt'), path('*_GenomePAM_raw.txt'), emit: genomePAM
    tuple val(meta), path('*PCV.txt'), emit: PCV

    script:
    """
    Rscript ${projectDir}/bin/genomePAM.R \
        ${identified_offtargets} ${projectDir}/resources/background_count/

    ## Filter by PAMlen, auto-detect position with highest enrichment
    best_pos=\$(awk -F'\t' 'NR>1 && \$6 == "${params.PAMlen}" {print \$4, \$14}' ${meta.id}_GenomePAM_raw.txt | sort -k2 -rn | head -1 | cut -d' ' -f1)
    awk -F'\t' -v pos="\$best_pos" '\$4 == pos && \$6 == "${params.PAMlen}" {print}' ${meta.id}_GenomePAM_raw.txt > ${meta.id}_PCV.txt
    """
}
