process genomePAM{
    tag{meta.id}

    conda params.r_conda

    input:
    tuple val(meta), path(identified_offtargets)

    output:
    path('*.html'), emit: genomePAM_html
    path('*.txt'), emit: genomePAM

    script:
    """
    ${params.r_conda}/bin/Rscript ${projectDir}/bin/genomePAM.R \
        ${identified_offtargets} ${projectDir}/resources/background_count/
    """
}