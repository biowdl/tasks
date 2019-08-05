Changelog
==========

<!--

Newest changes should be on top.

This document is user facing. Please word the changes in such a way
that users understand how the changes affect the new version.
-->

version 1.0.0-dev
---------------------------
+ GATK: Add CombineVariants task that allows, e.g., to merge VCFs from different callers.
+ Mutect2: Add GATK tasks related to variant filtering (LearnReadOrientationModel, MergeStats, GetPileupSummaries, CalculateContamination and FilterMutectCalls).
+ Mutect2: Add "--germline-resource" and "--f1r2-tar-gz" inputs, requiring an update to GATK 4.1.2.0. 
+ Mutect2: Add necessary missing index attribute for panel of normals.
+ VarDict: Add user definable flags (-M, -A, -Q, -d, -v, -f) to the paired VCF filtering script.
+ Cutadapt: If the output is a gzipped file, compress with level 1 (instead of default 6).
+ Cutadapt: Fix issues with read2output when using single-end reads.
+ Add feature type, idattr and additional attributes to htseq-count.
+ Added allow-contain option to bowtie.
+ Added a changelog to keep track of changes.
