#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include { FASTQC as FASTQC_PRE; FASTQC as FASTQC_POST } from './modules/nf-core/fastqc/main.nf'
include { MULTIQC as MULTIQC_PRE; MULTIQC as MULTIQC_POST } from './modules/nf-core/multiqc/main.nf'
include { trim_tag_umi } from  './modules/local/trim_tag_umi.nf'
include { trim_tag_umi_autodetect } from  './modules/local/trim_tag_umi_autodetect.nf'
include { align_identify } from  './modules/local/align_identify.nf'
include { annotate } from  './modules/local/annotate.nf'
include { visualize; consolidate_stat } from  './modules/local/visualize.nf'
include { svg_visualize } from './modules/local/svg_visualize.nf'
include { genomePAM } from './modules/local/genomePAM.nf'

// Create input channel from paired FASTQ files
ch_raw_short_reads = Channel.fromFilePairs(params.FQDIR + '/*_{R1,R2}*.fastq.gz', size: 2)
    .map { row ->
        def meta = [:]
        meta.id         = row[0].split(/[_\.]R1/)[0]
        meta.group      = 0
        meta.single_end = false
        return [ meta, row[1] ]
    }

// Reference genome: FASTA + BWA index files staged together
ch_genome_fasta = file(params.genome_fasta, checkIfExists: true)
ch_genome_index = Channel.fromPath("${params.genome_fasta}.{amb,ann,bwt,pac,sa,fai}", checkIfExists: true).collect()

workflow {
    // QC on raw sequencing reads
    FASTQC_PRE(ch_raw_short_reads)
    MULTIQC_PRE(FASTQC_PRE.out.zip.collect())

    // Trim adaptors + UMI extraction: cutadapt (known FIXSEQ) or BBDuk (auto-detect)
    if (params.FIXSEQ == "auto") {
        trim_tag_umi_autodetect(ch_raw_short_reads)
        ch_umitagged = trim_tag_umi_autodetect.out.umitagged_reads
        ch_num_reads = trim_tag_umi_autodetect.out.num_reads
    } else {
        trim_tag_umi(ch_raw_short_reads)
        ch_umitagged = trim_tag_umi.out.umitagged_reads
        ch_num_reads = trim_tag_umi.out.num_reads
    }

    // QC on trimmed + UMI-tagged reads
    FASTQC_POST(ch_umitagged)
    MULTIQC_POST(FASTQC_POST.out.zip.collect())

    // Align reads to reference genome + identify off targets
    align_identify(ch_umitagged, ch_genome_fasta, ch_genome_index)

    // Annotate off target sites (hg38 only) and create SVG visualizations
    if (params.GENOME == "hg38") {
        annotate(align_identify.out.identified_offtargets)
        svg_visualize(annotate.out.annotated_offtargets)
        genomePAM(annotate.out.annotated_offtargets)
        offtargets_for_viz = annotate.out.annotated_offtargets
    } else {
        svg_visualize(align_identify.out.identified_offtargets)
        offtargets_for_viz = align_identify.out.identified_offtargets
    }

    // Join channels for visualization:
    // [meta, offtargets] + [meta, umitagged]
    joined_ch = offtargets_for_viz
        .join(ch_umitagged)

    // Visualize the off target sites and sequence logos
    visualize(joined_ch)
    consolidate_stat(visualize.out.count_stat.collect())
}
