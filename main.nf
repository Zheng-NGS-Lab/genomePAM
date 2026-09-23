#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules

include { FASTQC as FASTQC_PRE; FASTQC as FASTQC_POST } from './modules/nf-core/fastqc/main.nf'  
include { MULTIQC as MULTIQC_PRE; MULTIQC as MULTIQC_POST } from './modules/nf-core/multiqc/main.nf'
include { trim_tag_umi } from  './modules/local/trim_tag_umi.nf'
include { consolidate } from  './modules/local/consolidate.nf'
include { align_identify } from  './modules/local/align_identify.nf'
include { annotate } from  './modules/local/annotate.nf'
include { visualize; consolidate_stat } from  './modules/local/visualize.nf'
include { svg_visualize } from './modules/local/svg_visualize.nf'
include { genomePAM } from './modules/local/genomePAM.nf'
include { chromatin_accessibility } from './modules/local/chromatin_accessibility.nf'

// log.info ""

workflow {
    // Report effective resource settings at startup
    log.info """
    genomePAM resource settings
      executor    : local
      max cpus    : ${params.max_cpus}
      max memory  : ${params.max_memory}
      bwa threads : ${params.BWATHREADS}
    """.stripIndent()

    // Validate that every param expected from the -params-file has a value
    // (list mirrors parameters.yml)
    def required_params = [
        FQDIR      : params.FQDIR,
        OUTDIR     : params.OUTDIR,
        Read1Tail  : params.Read1Tail,
        Read2Tail  : params.Read2Tail,
        pos1       : params.pos1,
        pos2       : params.pos2,
        posR2      : params.posR2,
        xNs        : params.xNs,
        FIXSEQ     : params.FIXSEQ,
        BWATHREADS : params.BWATHREADS,
        GENOME     : params.GENOME,
        AssaySpec  : params.AssaySpec,
    ]
    def empty_params = required_params.findAll { name, value ->
        value == null || value.toString().trim().isEmpty()
    }
    if (empty_params) {
        empty_params.each { name, value ->
            log.error "Param '${name}' has no value"
        }
        error "Found ${empty_params.size()} param(s) without a value in the -params-file -- fill in every field of the params file, then rerun"
    }

    // BWATHREADS / pos1 / pos2 / posR2 must be numbers (digits only)
    def numeric_params = [
        BWATHREADS : params.BWATHREADS,
        pos1       : params.pos1,
        pos2       : params.pos2,
        posR2      : params.posR2,
    ]
    def non_numeric_params = numeric_params.findAll { name, value ->
        !(value.toString().trim() ==~ /\d+/)
    }
    if (non_numeric_params) {
        non_numeric_params.each { name, value ->
            log.error "Param '${name}' is not a number: ${value}"
        }
        error "Found ${non_numeric_params.size()} non-numeric param(s) -- BWATHREADS, pos1, pos2 and posR2 must be numbers (digits only); fix the params file, then rerun"
    }

    // Warn (instead of aborting) when a genome other than hg38 is requested,
    // and state which analysis branches will run
    if (params.GENOME == 'hg38') {
        log.info "params.GENOME is hg38 -- annotate, svg_visualize and genomePAM will run"
    } else {
        log.warn "params.GENOME is '${params.GENOME}' (not hg38) -- annotate and genomePAM will skip, only svg_visualize will run"
    }

    if ((params.BWATHREADS as Integer) > (params.max_cpus as Integer)) {
        error "params.BWATHREADS (${params.BWATHREADS}) exceeds the executor cpu limit (${params.max_cpus}) -- align_identify tasks would be rejected; rerun with --BWATHREADS ${params.max_cpus} or lower"
    }

    // Validate required paths defined in nextflow.config before starting
    def required_paths = [
        genomePAM   : params.genomePAM,
        hg38        : params.hg38,
    ]
    def missing_paths = required_paths.findAll { name, path -> !file(path).exists() }
    if (missing_paths) {
        missing_paths.each { name, path ->
            log.error "Missing required path '${name}': ${path}"
        }
        error "Found ${missing_paths.size()} missing required path(s) defined in nextflow.config -- fix the params or create the files/directories, then rerun"
    }

    ch_raw_short_reads = Channel.fromFilePairs(params.FQDIR +'/*_{R1,R2}*.fastq.gz', size: 2)
    .map {
        row -> 
            def meta = [:]
            meta.id           = row[0].split(/[_\.]R1/)[0]
            meta.group        = 0
            meta.single_end   = false
            return [ meta, row[1] ]
    }

    // QC on raw sequencing reads
    FASTQC_PRE(ch_raw_short_reads)
    MULTIQC_PRE(FASTQC_PRE.out.zip.collect())

    // Trim i7 adaptors for read 1. Add UMI taggs to FASTQ
    trim_tag_umi(ch_raw_short_reads)

    // Consolidate reads from PCR duplicates. Record the number of reads
    consolidate(trim_tag_umi.out.umitagged_reads)

    // QC on trimmed + consolidated reads
    FASTQC_POST(consolidate.out.consolidated_reads)
    MULTIQC_POST(FASTQC_POST.out.zip.collect())

    // Align reads to human reference + identify off targets
    align_identify(consolidate.out.consolidated_reads)

    // Annotate off target sites only if aligned to hg38 and annotate into svg
    if (params.GENOME=="hg38") {
        annotate(align_identify.out.identified_offtargets)
        svg_visualize(annotate.out.annotated_offtargets)
        genomePAM(annotate.out.annotated_offtargets)
        // chromatin_accessibility(annotate.out.annotated_offtargets)
    } else {
        svg_visualize(align_identify.out.identified_offtargets)
    }

    // Join channels
    tmp_ch = align_identify.out.identified_offtargets.join(trim_tag_umi.out.umitagged_reads)
    // tmp_ch2 = tmp_ch.join(genomePAM.out.genomePAM)
    joined_ch = tmp_ch.join(consolidate.out.consolidated_reads)

    //Visualize the off target sites and sequence logos
    visualize(joined_ch)
    consolidate_stat(visualize.out.count_stat.collect())
}
