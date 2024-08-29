process chromatin_accessibility{
    tag {meta.id}

    conda params.r_conda
    errorStrategy 'ignore'
    input:
        tuple val(meta), path(identified_offtargets)

    output:
        path('*.pdf'), emit: chromatin_accessibility
    script:
    """
    echo "[LOG] Running chromatin_accessibility..."
    ${params.r_conda}/bin/Rscript ${projectDir}/bin/chromatin_accessibility.R \
        ${identified_offtargets}
    ## Move files to base directory
    mv output/*.pdf .
    rm output -r
    """
}