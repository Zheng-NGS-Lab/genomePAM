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

Channel.fromPath("${params.FQDIR}/**gz").set{fastq_gz_pair_ch}
ch_raw_short_reads = Channel.fromFilePairs(params.FQDIR +'/*_{R1,R2}*.fastq.gz', size: 2)
.map {
    row -> 
        def meta = [:]
		meta.id           = row[0].split(/[_\.]R1/)[0]
		meta.group        = 0
		meta.single_end   = false
		return [ meta, row[1] ]
}

workflow {
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
        chromatin_accessibility(annotate.out.annotated_offtargets)
    } else {
        svg_visualize(align_identify.out.identified_offtargets)
    }

    // Join channels
    tmp_ch = align_identify.out.identified_offtargets.join(trim_tag_umi.out.umitagged_reads)
    tmp_ch2 = tmp_ch.join(genomePAM.out.genomePAM)
    joined_ch = tmp_ch2.join(consolidate.out.consolidated_reads)

    //Visualize the off target sites and sequence logos
    visualize(joined_ch)
    consolidate_stat(visualize.out.count_stat.collect())
}
