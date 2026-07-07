process annotate{
    tag {meta.id}

    conda params.genomePAM

    input:
    tuple val(meta), path(identified_offtargets)

    output:
    tuple val(meta), path('*identifiedOfftargets.annotated.txt'), emit: annotated_offtargets

    beforeScript "mkdir identified"
    afterScript "rm -rf ./identified"

    script:
    """
    cp ${identified_offtargets} identified/
    # Annotate target sites with snpEff
    bash ${params.GUIDESEQDIR}/target_annotation_snpEff.sh \
        ${meta.id} ${params.GENOME} ${params.SNPEFFDIR}
    
    # Clean-up
    mv ./identified/*_identifiedOfftargets.annotated.txt .
    """
}
