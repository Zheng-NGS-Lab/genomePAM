process svg_visualize{
    tag {meta.id}

    input:
    tuple val(meta), path(offtargets)

    output:
    tuple val(meta), path('*.svg'), emit: offtarget_visual
    tuple val(meta), path('*.csv'), emit: offtarget_data

    script:
    """
    TargetSeq=\$(echo "${params.AssaySpec}" | tr '[:lower:]' '[:upper:]' | sed 's/_//g')

    # Determine PAM direction and extract PAM/spacer
    N="N"
    rest=\${TargetSeq#*\$N}
    N_idx=\$(( \${#TargetSeq} - \${#rest} - \${#N} ))
    if [[ \$N_idx -lt \$(( \${#TargetSeq} / 2 )) ]]; then
        PAM=\$(echo ${params.AssaySpec} | cut -d'_' -f1)
        TargetSeq_noPAM=\$(echo ${params.AssaySpec} | cut -d'_' -f2)
    else
        PAM=\$(echo ${params.AssaySpec} | cut -d'_' -f2)
        TargetSeq_noPAM=\$(echo ${params.AssaySpec} | cut -d'_' -f1)
    fi

    # Run guideseq.py visualize
    python3 ${params.GUIDESEQDIR}/guideseq.py visualize \
        --infile ${offtargets} \
        --PAM \${PAM^^} --target_seq \${TargetSeq^^} --outfolder . \
        --title ${meta.id}

    # Clean-up
    mv visualization/* .
    rmdir visualization
    """
}
