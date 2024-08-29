process visualize{
    tag {meta.id}

    conda params.r_conda
    errorStrategy 'ignore'
    input:
        tuple val(meta), path(identified_offtargets), path(umitagged_reads), path(consolidated_reads)

    output:
        tuple val(meta), path('*.pdf'), emit: sequence_logo
        path('.*_stats.csv'), emit: count_stat
        path('*.html'), emit: background_count
        
    script:
    """
    #####Get read counts
    umitagged_count=`grep '^@' ${umitagged_reads[0]} | wc -l`
    consolidated_count=`zgrep '^@' ${consolidated_reads[0]} | wc -l`

    #####Get TargetSeq and PAM from sample info CSV
    TargetSeq=\$(echo "${params.AssaySpec}"| tr '[:lower:]' '[:upper:]'|sed 's/_//g')
    # PAM: find the index of "N" and determine if PAM is before or after '_'
    N="N"
    rest=\${TargetSeq#*\$N}
    N_idx=\$(( \${#TargetSeq} - \${#rest} - \${#N} ))
    if [[ \$N_idx -lt \$(( \${#TargetSeq} / 2 )) ]]; then
        PAM=`echo ${params.AssaySpec} | cut -d'_' -f1`
        TargetSeq_noPAM=`echo ${params.AssaySpec} | cut -d'_' -f2`
    else
        PAM=`echo ${params.AssaySpec} | cut -d'_' -f2`
        TargetSeq_noPAM=`echo ${params.AssaySpec} | cut -d'_' -f1`
    fi

    # Run plot-PAM.R for visualization
    echo "[LOG] Running plot-PAM.R..."
    ${params.r_conda}/bin/Rscript ${projectDir}/bin/plot-PAM.R \
        ${identified_offtargets} \
        \$TargetSeq_noPAM \
        \$TargetSeq \
        ${meta.id} \$umitagged_count \$consolidated_count

    # Run plot-PAM_pos4.R to visualize the combinations
    echo "[LOG] Running plot-PAM_pos4.R..."
    run_date=\$(date "+%Y%m%d")
    
    ${params.r_conda}/bin/Rscript ${projectDir}/bin/plot-PAM_pos4.R \
        ${identified_offtargets} \$run_date
        
    echo "[LOG] Running counts script.R..."
    ${params.r_conda}/bin/Rscript ${projectDir}/bin/genomePAM.R \
        ${identified_offtargets} ${projectDir}/resources/background_count/
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
    echo "LibID,Umitagged_read_count,Unique_read_count,PM_n_reads,PM_n_sites,MM_n_reads,MM_n_sites" > allLib_stats.csv

    # Append content
    cat .*_stats.csv | grep -v '^LibID' | sort >> allLib_stats.csv
    """
}
