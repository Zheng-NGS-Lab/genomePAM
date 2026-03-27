# Changelog

All notable changes to the genomePAM pipeline are documented here.

## [1.2.0] - 2026-02-21

Major refactor replacing all custom UMI handling with industry-standard umi-tools.

### Changed

- **trim_tag_umi.nf**: Rewritten to use cutadapt (single-pass paired-end adapter trimming) and umi-tools extract (UMI extraction with barcode whitelist filtering). Replaces three sequential BBDuk calls and custom umitag.py. cutadapt was chosen over BBDuk for the known-FIXSEQ path because it handles paired-end trimming and coordinated length filtering in a single pass, with structured logs consumed by MultiQC.
- **align_identify.nf**: Added `umi-tools dedup --paired --method directional` after BWA alignment. Deduplication moved from pre-alignment (FASTQ-level molecular ID grouping) to post-alignment (genomic coordinate + UMI with network-based clustering). This eliminates false duplicate calls from reads with the same UMI mapping to different genomic locations.
- **main.nf**: Removed consolidate process. Added if/else routing: `params.FIXSEQ == "auto"` selects BBDuk autodetect path; explicit FIXSEQ selects cutadapt path. Downstream channels unified via `ch_umitagged`.
- **identifyOfftargetSites.py**: Updated parseReadName() to parse umi-tools `@READ_UMI` header format.
- **visualize.nf**: Removed consolidated_reads input channel (no longer produced).
- **plot-PAM.R**: Removed consolidated_count parameter from plots.
- **Dockerfile**: Added cutadapt and umi_tools to mamba install.
- **Nextflow script escaping**: Normalized all bash line continuations to `\\` (proper Nextflow GString idiom) instead of mixing single `\` (Groovy line continuation).

### Added

- **trim_tag_umi_autodetect.nf**: New process for FIXSEQ auto-detection. Preserves the BBDuk-based sequential trimming approach since the barcode must be discovered from R1 data before computing its reverse complement for R2 trimming. Includes reformat.sh paired-end length filter before umi-tools extract to prevent crashes on short reads.
- **tests/integration/**: Integration test suite covering off-target identification, FASTQ reader validation, read name parsing, alignment sequences, and real data structural properties.
- **tests/test_fastq_reader.py**: Unit tests validating fq() reader against real downsampled data.
- **tests/data/fastq/**: Downsampled 10K-read FASTQ fixtures (cas9 SRR33421097, fncas12a GNPAM327, 3.5 MB total). Committed directly to git rather than downloaded in CI to avoid SRA network dependencies and keep tests self-contained.
- **tests/compare_trim_approaches.py**: Comparison utility validating cutadapt vs BBDuk produce equivalent output.
- **conf/modules.config**: Added resource configuration for trim_tag_umi_autodetect process.

### Removed

- **modules/guideseq/umi/umitag.py**: Custom FASTQ reader/writer with umitag_inline(), build_molecular_id(), parse_whitelist(). Replaced by umi-tools extract which is the standard tool for this operation and handles edge cases the custom code did not (crashes on reads shorter than 18bp barcode+UMI pattern).
- **modules/guideseq/umi/consolidate.py**: Pre-alignment FASTQ deduplication by molecular ID. Replaced by post-alignment umi-tools dedup which uses both genomic coordinates and UMI for more accurate duplicate detection.
- **modules/local/consolidate.nf**: Nextflow process wrapper for consolidate.py.
- **tests/test_umitag.py**, **tests/test_umitag_inline.py**, **tests/test_consolidate.py**: Tests for deleted modules.

## [1.1.0] - 2025-08-15

UMI handling improvements and code cleanup.

### Changed

- **umitag.py**: Replaced fragile bash/awk UMI extraction with Python inline processing mode (umitag_inline). Produces umi-tools compatible headers (`@READID_UMI 1:N:0:0`).
- **identifyOfftargetSites.py**: Auto-detects PCV (position-count-value) column from data instead of using PAMpos as a literal filter. Enables correct operation with both 3' and 5' PAM nucleases.
- **guideseq.py**: Removed dead code, fixed deprecated file modes (`'rU'` to `'r'`), closed resource leaks.

## [1.0.0] - 2025-06-01

Initial release. Containerized fork with end-to-end pipeline fixes.

### Added

- **Dockerfile**: Single container with all dependencies (BWA, samtools, BBMap, snpEff, Python 3.10, R 4.3). Replaces Conda-based setup.
- **run.sh**: Automated setup script (Docker build, data download, reference indexing, pipeline execution).
- **FIXSEQ auto-detection**: Detects conserved barcode sequence from data at R1 positions pos1-pos2. Replaces hardcoded SpCas9 barcode value.
- **CI**: pytest workflow on merge to main.

### Fixed

- **SRA FASTQ header support**: Fixed parsing of 3-field SRA-format headers (`@SRR... length=N`) that lacked instrument/flowcell fields expected by the original BBDuk-based trimming.
- **Annotation column alignment**: Added empty trailing tab columns for unannotated rows so all rows have the same column count, preventing `read.table` parse failures in R.
- **R script numeric coercion**: Fixed `cumsum()`, `rep()`, and `sum()` on character data caused by `colClasses=c("character")` in `read.table`. Added `as.numeric()` coercion and `fill=TRUE` for defensive TSV parsing.
- **Vendored UMI package**: Restored `umi` Python package (from `aryeelab/umi`) into `modules/guideseq/umi/`, fixing broken imports in guideseq.py.
