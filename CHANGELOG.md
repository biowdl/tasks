Changelog
==========

<!--

Newest changes should be on top.

This document is user facing. Please word the changes in such a way
that users understand how the changes affect the new version.
-->

version 1.0.0-dev
---------------------------
+ GATK: SplitNCigarReads, BaseRecalibration and ApplyBQSR do no longer need regions files as required inputs.
+ Cutadapt: If the output is a gzipped file, compress with level 1 (instead of default 6).
+ Cutadapt: Fix issues with read2output when using single-end reads.
+ Add feature type, idattr and additional attributes to htseq-count.
+ Added allow-contain option to bowtie.
+ Added a changelog to keep track of changes.
