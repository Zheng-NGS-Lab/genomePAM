process genomePAM{
    tag{meta.id}

    conda params.r_conda
    errorStrategy 'ignore'
    input:
    tuple val(meta), path(identified_offtargets)

    output:
    path('*.html'), emit: genomePAM_html
    path('*.txt'), emit: genomePAM
    path('*PCV.txt'), emit: PCV
    
    script:
    """
    echo "[LOG] Running counts script.R..."
    ${params.r_conda}/bin/Rscript ${projectDir}/bin/genomePAM.R \
        ${identified_offtargets} ${projectDir}/resources/background_count/

    ## Filter by PAMlen and Position
    awk -F'\t' '\$4 == "${params.PAMpos}" && \$6 == "${params.PAMlen}" {print}' ${meta.id}_GenomePAM_raw.txt > ${meta.id}_PCV.txt

    """
}
