Changelog
==========

<!--

Newest changes should be on top.

This document is user facing. Please word the changes in such a way
that users understand how the changes affect the new version.
-->

version 1.1.0-dev
---------------------------
+ TALON: Update average memory needs for main TALON process

version 1.0.0
---------------------------
+ Common: Add "SampleConfigToSampleReadgroupLists" task
+ MultiQC: the "interactive" input is now set to true by default
+ Removed deprecated tasks:
  + bioconda.installPrefix
  + mergecounts.MergeCounts
+ GATK.BaseRecalibrator: "knownIndelsSitesVCFs" and "knownIndelsSitesVCFIndexes" are no longer optional, but now have a default of "[]"
+ Removed BWA index task
+ Removed unused "picardJar" input from bwa.wdl
+ All inputs to bedtools Sort are now reflected in the generated command
+ TranscriptClean: Update TranscriptClean container to version 1.0.8
+ Removed "pipefail" from command sections TALON and TranscriptClean
+ Add WDL task for Minimap2
+ Add WDL task for TALON
+ Add WDL task for TranscriptClean
+ Fastqsplitter: fix mkdir command to work with biocontainer's busybox mkdir
+ Cutadapt: simplify interface
+ Bigger memory multiplier in mutect to take in account bigger vmem usage
+ Cutadapt: Remove default adapter
+ Fastqsplitter: use version 1.1.
+ Picard: Use version 2.20.5 of the biocontainer as this includes the R dependency
+ Common: Update dockerTag to dockerImage.
+ GATK: Add CombineVariants task that allows, e.g., to merge VCFs from different callers.
+ Mutect2: Add GATK tasks related to variant filtering (LearnReadOrientationModel, MergeStats, GetPileupSummaries, CalculateContamination and FilterMutectCalls).
+ Mutect2: Add "--germline-resource" and "--f1r2-tar-gz" inputs, requiring an update to GATK 4.1.2.0. 
+ Mutect2: Add necessary missing index attribute for panel of normals.
+ MultiQC: Add memory variable to multiqc task.
+ GATK: SplitNCigarReads, BaseRecalibration and ApplyBQSR do no longer need regions files as required inputs.
+ VarDict: Add user definable flags (-M, -A, -Q, -d, -v, -f) to the paired VCF filtering script.
+ Cutadapt: If the output is a gzipped file, compress with level 1 (instead of default 6).
+ Cutadapt: Fix issues with read2output when using single-end reads.
+ Add feature type, idattr and additional attributes to htseq-count.
+ Added allow-contain option to bowtie.
+ Added a changelog to keep track of changes.
+ Added sortByName task in samtools to support more memory efficient execution of HTSeqCount.
+ Removed the bam index from HTSeqCount's inputs.
