# WGS-Somatic-Variant-Pipeline

## Workflow

* 1.WGSseq.map.R - Format and concatenate fastq files
* 2.WGSseq.map.R - FastqToSam: Converting Fastq to uBAM
* 3.1.WGSseq.map.R - SamToFastqAndBwaMem: Read unmapped BAM, convert on-the-fly to FASTQ and stream to BWA MEM for alignment
* 3.2.WGSseq.map.R - MergeBamAlignment: Merge original input uBAM file with BWA-aligned BAM file
* 3.3.WGSseq.map.R - MarkDuplicates: Mark duplicate reads to avoid counting non-independent observations
* 3.4.WGSseq.map.R - SortAndFixTags: Sort BAM file by coordinate order and fix tag values for NM and UQ
* 3.5.WGSseq.map.py - create sequence_grouping.txt and sequence_grouping_with_unmapped.txt
* 3.6.WGSseq.map.R - BaseRecalibrator: Perform Base Quality Score Recalibration (BQSR) on the sorted BAM in parallel
* 3.7.WGSseq.map.R - GatherBQSRReports: Merge/gather the recalibration reports resulting from by-interval recalibration
* 3.8.WGSseq.map.R - ApplyBQSR: Apply Base Quality Score Recalibration (BQSR) model
* 3.9.WGSseq.map.R - GatherBamFiles: Merge the recalibrated BAM files resulting from by-interval recalibration
* 4.0.SplitIntervals.R - splits the master interval file into interval files for scattering resulting files contain equal number of bases.
* 4.1.CreateSomaticPON.R - Mutect2: creates vcf files for each normal sample
* 4.2.CreateSomaticPON.R - CreateSomaticPanelOfNormals: Create PON databse and this script combines the normal calls using CreateSomaticPanelOfNormals