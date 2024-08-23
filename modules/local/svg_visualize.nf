process svg_visualize{
    tag {meta.id}

    conda params.guide_seq_conda

    input:
        tuple val(meta), path(offtargets)

    output:
    tuple val(meta), path('*.svg'), emit: offtarget_visual

    script:
    """
    # Add GUIDEseq env bin to PATH
    PATH=${params.guide_seq_conda}/bin:\$PATH
    
    #####Get TargetSeq and PAM from sample info CSV
    TargetSeq=\$(echo "${params.AssaySpec}"| tr '[:lower:]' '[:upper:]'|sed 's/_//g')

    # TargetSeq
    # PAM: find the index of "N" and determine if PAM is before or after '_'
    N="N"
    rest=\${TargetSeq#*\$N}
    N_idx=\$(( \${#TargetSeq} - \${#rest} - \${#N} ))
    if [[ \$N_idx -lt \$(( \${#TargetSeq} / 2 )) ]]; then
        PAM=`echo ${params.AssaySpec} | cut -d'_' -f1`
        TargetSeq_noPAM=`echo \$AssaySpec | cut -d'_' -f2`
    else
        PAM=`echo ${params.AssaySpec} | cut -d'_' -f2`
        TargetSeq_noPAM=`echo \$AssaySpec | cut -d'_' -f1`
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