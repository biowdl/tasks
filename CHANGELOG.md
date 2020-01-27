Changelog
==========

<!--

Newest changes should be on top.

This document is user facing. Please word the changes in such a way
that users understand how the changes affect the new version.
-->

version 2.2.0-dev
---------------------------
+ Increase default memory of BWA mem to 32G (was 16G).
+ Add common.TextToFile task.
+ Add bedtools.Intersect.
+ Add `-o pipefail` to bedtools.MergeBedFiles to prevent errors in BED files 
  from going unnoticed.
+ Centrifuge: Fix -1/-U options for single end data.
+ Add bedtools.Complement, bedtools.Merge, and add a task to combine multiple
  bed files called bedtools.MergeBedFiles. This task combines bedtools merge 
  and sort.
+ Change `g` parameter on bedtools.Sort to `genome`.
+ Add `ploidity` and `excludeIntervalList` to gatk.HaplotypeCallerGvcf.
+ Update centrifuge tasks.
+ Removed unused "cores" inputs from transcriptclean tasks.
+ Removed unused "cores" inputs from talon tasks.
+ Removed unused "threads" input from ModifyStrelka.
+ Removed the "installDir" inputs from the somaticseq tasks.
+ Removed the "installDir" input from CombineVariants.
+ Removed the "extraArgs" input from FilterMutectCalls.
+ Removed unused "verbose" and "quiet" inputs from multiqc.
+ Added parameter_meta sections to a variety of tasks.
+ Picard's BedToIntervalList outputPath input is now optional (with a default of "regions.interval_list")
+ TALON: Fix SQLite error concerning database/disk space being full.
+ Update htseq to default image version 0.11.2
+ Update biowdl-input-converter in common.wdl to version 0.2.1.
+ Update TALON section to now include the new annotation file output, and add config file creation to the TALON task.
+ Removed unused inputs (trimPrimer and format) for cutadapt.
+ Various minor command tweaks to increase stability.
+ Fixed unused inputs in bedtools sort (inputs are now used).
+ Added miniwdl check to linting.
+ Update TALON default image to version 4.4.1.

version 2.1.0
---------------------------
+ Make intervals optional for GATK CombineGVCFs.
+ Updated biowdl-input-converter version.
+ GATK CombineGVCFs memory was tripled to prevent it from using a lot of CPU in
  Garbage Collection mode.
+ Updated parameter_meta sections for Minimap2 and TranscriptClean to wdl-aid format.
+ Updated cores variable for TALON, the default is now 4.
+ Updated TALON to version 4.4.
+ Added parameter_meta sections to the following tools:
    + htseq
    + cutadapt
    + collect-columns
    + stringtie
    + fastqc
+ Updated star default image to 2.7.3a.
+ Hisat2 now indexes the resulting BAM file.
+ Samtools index now also works without setting a path for the output.
+ Bugfix: Biowdl-input-converter now makes sure the output directory exists.

version 2.0.0
---------------------------
+ TranscriptClean: Update TranscriptClean to version 2.0.2.
+ Memory runtime attributes are now Strings indicating total memory, as opposed to Ints indicating memory per core.
+ Memory inputs for most tasks are now Strings, remaining Int memory inputs are renamed to "memoryGb".
+ Use the biowdl-input-converter container for JsonToYaml, to reduce the amount of containers needed.
+ Add biowdl-input-converter and remove SampleConfigToSampleReadgroupLists which it replaces.
+ GATK.GenotypeGVCFs: Increased memoryMultiplier from 2.0 to 3.0 .
+ Minimap2: Add -k option to minimap2 mapping.
+ Added bwakit task.
+ Minimap2: Add the option for --MD tag.
+ TALON: Update average memory needs for main TALON process.

version 1.0.0
---------------------------
+ Common: Add "SampleConfigToSampleReadgroupLists" task.
+ MultiQC: the "interactive" input is now set to true by default.
+ Removed deprecated tasks:
  + bioconda.installPrefix
  + mergecounts.MergeCounts
+ GATK.BaseRecalibrator: "knownIndelsSitesVCFs" and "knownIndelsSitesVCFIndexes" are no longer optional, but now have a default of "[]".
+ Removed BWA index task.
+ Removed unused "picardJar" input from bwa.wdl.
+ All inputs to bedtools Sort are now reflected in the generated command.
+ TranscriptClean: Update TranscriptClean container to version 1.0.8.
+ Removed "pipefail" from command sections TALON and TranscriptClean.
+ Add WDL task for Minimap2.
+ Add WDL task for TALON.
+ Add WDL task for TranscriptClean.
+ Fastqsplitter: fix mkdir command to work with biocontainer's busybox mkdir.
+ Cutadapt: simplify interface.
+ Bigger memory multiplier in mutect to take in account bigger vmem usage.
+ Cutadapt: Remove default adapter.
+ Fastqsplitter: use version 1.1.
+ Picard: Use version 2.20.5 of the biocontainer as this includes the R dependency.
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
