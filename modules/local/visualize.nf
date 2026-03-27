process visualize{
    tag {meta.id}

    errorStrategy 'ignore'

    input:
    tuple val(meta), path(identified_offtargets), path(umitagged_reads)

    output:
    tuple val(meta), path('*.pdf'), emit: sequence_logo
    path('.*_stats.csv'), emit: count_stat

    script:
    """
    #####Get read counts
    umitagged_count=\$(zcat ${umitagged_reads[0]} | grep '^@' | wc -l)

    #####Get TargetSeq and PAM from AssaySpec
    TargetSeq=\$(echo "${params.AssaySpec}" | tr '[:lower:]' '[:upper:]' | sed 's/_//g')
    # PAM: find the index of "N" and determine if PAM is before or after '_'
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

    # Run plot-PAM.R for sequence logo visualization
    Rscript ${projectDir}/bin/plot-PAM.R \
        ${identified_offtargets} \
        \$TargetSeq_noPAM \
        \$TargetSeq \
        ${meta.id} \$umitagged_count

    # Run plot-PAM_pos4.R for 4-position heatmaps
    run_date=\$(date "+%Y%m%d")
    Rscript ${projectDir}/bin/plot-PAM_pos4.R \
        ${identified_offtargets} \$run_date
    """
}

process consolidate_stat {
    input:
    path count_stat_csvs

    output:
    path "allLib_stats.csv", emit: consolidated_stat

    script:
    """
    # Write header
    echo "LibID,Umitagged_read_count,PM_n_reads,PM_n_sites,MM_n_reads,MM_n_sites" > allLib_stats.csv

    # Append content
    cat .*_stats.csv | grep -v '^LibID' | sort >> allLib_stats.csv
    """
}
