Changelog
==========

<!--
Newest changes should be on top.

This document is user facing. Please word the changes in such a way
that users understand how the changes affect the new version.
-->
version 5.1.0-dev
---------------------------
+ Purple's `somaticRainfallPlot` output is now optional and included in
  the `plots` output as well.
+ Bedtools coverage's timeMinutes now defaults to `320`.
+ Gridss' runtime attribute defaults were changed to:
  + jvmHeapSizeGb: `64`
  + nonJvmMemoryGb: `10`
  + threads: `12`
+ Virusbreakend's runtime attribute defaults were changed to:
  + threads: `12`
  + timeMinutes: `320`
+ Cobalt's timeMinutes now defaults to `480`.
+ Orange's timeMinutes now defaults to 10.
+ Sage's runtime attributes were changed to:
  + threads: `32`
  + javaXmx: `"16G"`
  + memory: `"20G"`
  + timeMinutes: `720`
+ Sambamba's runtimeMinutes nor defaults to `320`.
+ Added a task for CupGenerateReport.
+ Updated Cuppa to version 1.6.
+ Added a task for Gripss.
+ Fixed the HealthChecker task's determination of the `succeeded` output
  value.
+ Updated Linx to version 1.18.
+ Added a task for LinxVisualization.
+ Added a task for HMFtools Orange.
+ Added a task for HMFtools Pave.
+ Updated Purple to version 3.2.
+ Added plot and table outputs of Sage to task outputs.
+ Updated virus-interpreter to version 1.2.
+ Updated Peach to version 1.5.
+ Added a task to add SVTYPE annotations to GRIDSS results
  (`AnnotateSvTypes`).
+ The GRIDSS task will now run tabix separately if GRIDSS doesn't
  produce a vcf index.
+ Add a script to subtract UMI's from the read name and add them as
  a BAM tag for each BAM record. The script is in umi.BamReadNameToUmiTag.
+ Add fgbio.AnnotateBamWithUmis.
+ Add picard.UmiAwareMarkDuplicatesWithMateCigar.
+ Added a task for SnpEff.
+ Adjusted runtime settings for sambamba Markdup.
+ Added a task for sambamba Flagstat.
+ Added a task for Picard CollectWgsMetrics.
+ Added a task for Peach.
+ Added tasks for HMFtools:
  + Amber
  + Cobalt
  + Cuppa
  + CuppaChart
  + GripssApplicationKt
  + GripssHardFilterApplicationKt
  + HealthChecker
  + Linx
  + Protect
  + Purple
  + Sage
  + VirusInterpreter
+ Added a task for VirusBreakend.
+ Added a task for GridssAnnotateVcfRepeatmasker.
+ Bumped GRIDSS version to 2.12.2.
+ Adjusted GRIDSS runtime settings.
+ Added optional inputs to GRIDSS:
  + blacklistBed
  + gridssProperties
+ Added a task for GRIDSS AnnotateInsertedSequence.
+ Added a task for ExtractSigPredictHRD.
+ Added a task for DeconstructSigs.
+ Added option useSoftclippingForSupplementary (default false) to
  BWA mem.
+ Adjusted BWA mem runtime settings.
+ Added a task for bedtools coverage.
+ Added a task for bcftools filter.
+ Adjusted runtime settings for bcftools annotate.
+ Added optional inputs to bcftools annotate:
  + inputFileIndex
  + annsFileIndex
+ Update parameter_meta for macs2
+ Add sample position in array task.

version 5.0.2
---------------------------
+ bumped ScatterRegions container to 1.0.0

version 5.0.1
---------------------------
+ Smoove: enable genotyping
+ add runtime memory to number of tasks.

version 5.0.0
---------------------------
+ Update CPAT to version 3.0.4.
  + Changed the `outFilePath` input to `outputPrefix`.
+ Survivor: Change integer to string literal in boolean parameters.
+ Samtools: Add mkdir line to `Fastq` task.
+ Add new parameters from CCS version 6.0.0 and add two new outputs:
  `ccs_report.txt` & `zmw_metrics.json.gz`.
+ Change CutAdapt memory to `5G`.
+ Increase multiqc base time from 5 to 10.
+ Update biowdl-input-converter to version 0.3.
+ Update minimap2 to version 2.20.
+ Update lima to version 2.2.0.
+ Update ccs to version 6.0.0.
+ Update bam2fastx to version 1.3.1.
+ Add memory values to GffCompare, GffRead and CPAT.
+ GffCompare: Make the `referenceAnnotation` input optional.
+ Stringtie: Add the `minimumCoverage` input.
+ UMI-tools: Update default dockerImage to use umitools v1.1.1 with correct
             samtools version (1.10).
+ UMI-tools: Re-introduce samtools indexing.
+ UMI-tools: Update default dockerImage to use umitools v1.1.1.
+ UMI-tools dedup: Add tempdir.
+ Bcftools view: Add options for filtering (include, exclude, excludeUncalled).
+ Duphold: Add `duphold.wdl`.
+ Add new wdl file prepareShiny.wdl for creating input files for shiny app.
+ mergePacBio: Rename `mergedReport` to `outputPathMergedReport`.
+ Lima: Fix copy commands.
+ Fixed the `size` call in the default for gffread's timeMinutes, to retrieve
  GBs instead of bytes.
+ Update stringtie to version 1.3.6.
+ Update Lima to version 2.0.0.
+ Update IsoSeq3 to version 3.4.0.
+ Update samtools to version 1.11.
+ Update Picard to version 2.23.8.
+ Update NanoPlot to version 1.32.1.
+ Update MultiQC to version 1.9.
+ ~Update StringTie to version 2.1.4.~
+ Complete `parameter_meta` for tasks missing the outputs.
+ DeepVariant: Add an optional input for the gvcf index.
+ Samtools: `Sort` task now has `threads` in runtime instead of `1`.
+ Picard: Add parameter_meta to `SortSam`.
+ pbmm2: Add parameter_meta for `sample`.
+ Centrifuge: Rename output in task `KReport` to `KrakenReport` to resolve
  name collision with task name.
+ Bwa & bwa-mem2: Add parameter_meta for `outputHla`.
+ Multiqc: Removed WDL_AID excludes of "finished" & "dependencies" inputs.
+ Bam2fastx: Add localisation of input files to Bam2Fasta task.
+ Lima: `cores` input has been renamed to `threads` to match tool naming.
+ isoseq3: `cores` input has been renamed to `threads` to match tool naming.
+ CCS: `cores` input has been renamed to `threads` to match tool naming.
+ Add PacBio preprocessing specific tasks `mergePacBio` & `ccsChunks`.
+ CCS: Update CCS to version 5.
+ deepvariant: Add task for DeepVariant.
+ gatk: Make intervals optional for GenotypeGVCFs.
+ isoseq3: Add required bam index input to isoseq3.
+ pbbam: Add task for indexing PacBio bam files.
+ picard: Add CollectHsMetrics and CollectVariantCallingMetrics.
+ Samtools: Add `threads` to parameter meta for Merge task.
+ bcftools: add tmpDir input to specify temporary directory when sorting.
+ bcftools: remove outputType and implement indexing based on output
  file extension.
+ NanoPack: Add parameter_meta to NanoPlot task.
+ Centrifuge: Remove metrics file from classification (which causes the
  summary report to be empty).
  https://github.com/DaehwanKimLab/centrifuge/issues/83
+ Add NanoPlot and NanoQC tasks.
+ Centrifuge: Add `timeMinutes` to `Classify` task and remove unnecessary
  downloading tasks (alternative is refseqtools).
+ collect-columns: updated docker image to version 1.0.0 and added the
  `sumOnDuplicateId` input (defaults to false).
+ survivor: replace integer boolean type to logical true or false value.
+ vt: Add option to ignore masked reference.
+ bcftools: add sorting and annotation.
+ Bam2fastx: Input bam and index are now arrays.
+ Lima: Remove globs from outputs.
+ Updated task gridss.wdl: add --jvmheap parameter.
+ A bwa-mem2 task was created with the same interface (including usePostalt)
  as the bwa mem task.
+ bwa mem and bwa kit are now one task. The usePostalt boolean can be used to
  switch the postalt script on and off.
+ Added a task for GRIDSS.
+ Add wdl file for pacbio's bam2fastx tool.

version 4.0.0
---------------------------
+ Picard MergeVcf now uses compression level 1 by default.
+ bwa mem, bwa mem+kit and hisat2 have their samtools sort threads tweaked. The
  number of threads is now related to the number of threads on the aligner.
  Using more threads reduces the chance of the samtools sort pipe getting
  blocked if it's full.
+ Renamed a few inputs in centrifuge.wdl, isoseq3.wdl, talon.wdl,
  transcriptclean.wdl to be more descriptive.
+ Renamed outputs of tasks used in the TALON-WDL, PacBio-subreads-processing &
  sequence-classification pipelines.
+ Reworked bcf2vcf task into bcftools view task.
+ Removed the redundant format flag from the htseq interface. This is
  autodetected in newer versions of htseq.
+ Update docker images for samtools, bcftools, picard, GATK, cutadapt, htseq
  and chunked-scatter.
+ Default docker images for bwa, bwakit and hisat2 updated to include samtools
  1.10.
+ Alignment tasks (STAR, Hisat2, BWA) now produce BAM files at level 1
  compression.
+ Hisat2 task has added controls for samtools.
+ Alignment tasks no longer produce BAM indexes as these are not needed
  by the markduplicates step.
+ Picard Markduplicates now uses 7G of RAM just like in GATK's best practice
  example pipeline.
+ Picard SortSam added as a task.
+ Md5 files are no longer created by default on Picard tasks that generate
  BAM files.
+ Changed PicardMarkduplicates to use COMPRESSION_LEVEL=1 by default with
  the htsjdk deflater.
  This makes the task finish in 32% less time at the cost of a 8% larger BAM
  file.
+ Added sambamba markdup and sambamba sort. NOTE: samtools sort is more
  efficient and is recommended.
+ Correctly represent samtools inconsistent use of the threads flag.
  Sometimes it means 'threads' sometimes it means 'additional threads'.
  BioWDL tasks now use only threads. The `threads - 1` conversion is
  applied where necessary for samtools tools that use additional threads.
+ Updated BWA MEM  and BWA KIT tasks to use samtools sort version 1.10 for
  sorting the BAM file.
+ Updated memory requirements on bcftools Stats, bwa mem, bwakit, GATK
  ApplyBQSR, GATK BaseRecalibrator, GATK GatherBqsrReports, Gatk
  HaplotypeCaller, Picard CollectMultipleMetrics, Picard GatherBamFiles,
  samtools Flagstat, samtools sort and bcftools stats.
+ TALON: Update `FilterTalonTranscripts` to new version, which removes the
  pairingsFile and replaces this with datasetsFile.
+ TALON: Add `GetSpliceJunctions` & `LabelReads` tasks.
+ TALON: Update to version 5.0.
+ Add tasks for pbmm2, the PacBio wrapper for minimap2.
+ Update the image for chunked-scatter and make use of new features from 0.2.0.
+ Tuned resource requirements for GATK VariantEval, MultiQC, Picard metrics and
  STAR.
+ Added a new task for [scatter-regions](https://github.com/biowdl/chunked-scatter)
  that replaces biopet-scatterregions.
+ The FastQC task now talks to the Java directly instead of using the included
  Perl wrapper for FastQC. This has the advantage that memory and threads can
  be set independently. A rather high maximum heap size of 1750MB (Xmx1750M)
  was set, as OOM errors occurred frequently on some fastqs.
+ STAR: Add options regarding alignment score (regarding read length as well)
  for tweaking when processing rRNA depleted samples.
+ TALON: Update `minimumIdentity` to correct type (float, was integer)
  & set new default according to developers (0.8, was 0).
+ Added GATK VariantEval task.
+ Added a log output for STAR.
+ Added report output to Hisat2.
+ Added output with all reports to gffcompare.
+ Change MultiQC inputs. It now accepts an array of reports files. It does not
  need access to a folder with the reports anymore. MultiQC can now be used
  as a normal WDL task without hacks.
+ Picard: Make all outputs in `CollectMultipleMetrics` optional. This will
  make sure the task will not fail if one of the metrics is set to false.
+ The struct `BowtieIndex` was removed, as it has become obsolete.
+ The task `ReorderGlobbedScatters` was removed, as it has become obsolete.
+ Adjusted the memory settings of many tools, especially java tools.
  They should now more accurately represent actual memory usage (as
  opposed to virtual memory).
+ Added `-XX:ParallelGCThreads=1` to the java options of java tasks.
+ Added `timeMinutes` input to many tasks, this indicates a maximum
  number of minutes that the job will run. The associated runtime
  attribute is `time_minutes` which can be used to inform
  a scheduler (eg. slurm) of the run time of the job.
+ Added STAR GenomeGenerate task.
+ GATK.HaplotypeCaller: Add `--dont-use-soft-clipped-bases` and
  `--standard-min-confidence-threshold-for-calling` options. These are
  required for RNA seq variant calling according to GATK best practices.
+ Samtools: Fix quotations in sort command.
+ Samtools SortByName is now called Sort.
+ Generalize sort task to now also sort by position, instead of just read name.
+ Add CreateSequenceDictionary task to picard.
+ Add faidx task to samtools.
+ Isoseq3: Remove dirname command from output folder creation step.
+ Isoseq3: Requires more memory by default, is now 2G.
+ Isoseq3: Remove cp commands and other bash magic, file naming is now
  solved by pipeline.
+ Lima: Replace mv command with cp.
+ Add WDL task for smoove (lumpy) sv-caller.

version 3.1.0
---------------------------
+ Default threads for BWA in bwa.Kit task: 4. Samtools sort in the
  same task: 1. Output BAM compression level to 1.
+ Lima: Add missing output to parameter_meta.
+ Lima: Remove outputPrefix variable from output section.
+ Isoseq3: Make sure stderr log file from Refine is unique and not overwritten.
+ Isoseq3: Add workaround in Refine for glob command not locating files
  in output directory.
+ Isoseq3: Fix --min-polya-length argument syntax.
+ Lima: Add workaround for glob command not locating files in output directory.
+ CCS: Add missing backslash.
+ Cutadapt now explicitly calls the `--compression-level` flag with compression
  level 1 to prevent cutadapt from using very high gzip compression level 6
  that uses 400% more cpu time.
+ Update default docker image for cutadapt and fastqc.
+ Default number of cores for cutadapt and bwamem to 4 cores.

version 3.0.0
---------------------------
+ Add optional input umiSeparator in umi-tools dedup task.
+ Update command section syntax Minimap2, Talon, TranscriptClean and Centrifuge.
+ Add CCS workflow WDL files (ccs.wdl, lima.wdl, isoseq3.wdl).
+ Update TALON version to 4.4.2.
+ The statsPrefix input for umitools dedup is now optional.
+ Allow setting the `--emit-ref-confidence` flag for HaplotypeCaller.
+ Add `--output-mode` flag to HaplotypeCaller.
+ Added rtg.Format and rtg.VcfEval tasks.
+ Added gatk.SelectVariants and gatk.VariantFiltration tasks.
+ Fixed a bug where the output directory was not created for bwa.Kit.
+ Add vt task for variants normalization and decomposition.
+ Update WDL task Picard (Add task RenameSample).
+ Update WDL task Samtools (Add task FilterShortReadsBam).
+ Add WDL task for BCFtools (bcf to vcf).
+ Add WDL task for SURVIVOR (merge).
+ Update WDL task Manta (Add germline SV calling).
+ Add WDL task for Delly.
+ Add WDL task for Clever (and Mate-Clever).
+ Add proper copyright headers to all WDL files. So the free software license
  is clear to end users who wish to adapt and modify.
+ Add pedigree input for HaplotypeCaller and GenotypeGVCFs.
+ Combined biopet.ScatterRegions and biopet.ReorderedGlobbedScatters into one.
  biopet.ScatterRegions now always returns correctly ordered scatters.
+ Add tasks for umi-tools dedup and extract.
+ Add `GenomicsDBImport` task for GATK.
+ Add `annotationGroups` input to `GenotypeGVCFs` to allow setting multiple
  annotation groups. The `StandardAnnotation` group is still used as default.
+ GenotypeGVCFs, only allow one input GVCF file, as the tool also only allows
  one input file.
+ Rename HaplotypeCallerGVCF to HaplotypeCaller. Add `gvcf` option to set
  whether output should be a GVCF.
+ Centrifuge: Add Krona task specific to Centrifuge.
+ Centrifuge: Fix Centrifuge tests, where sometimes the index files could
  still not be located.
+ Update parameter_meta for TALON, Centrifuge and Minimap2.
+ Centrifuge: Fix issue where Centrifuge Inspect did not get the correct
  index files location.
+ Add `minimumContigLength` input to PlotDenoisedCopyRatios
  and PlotModeledSegments.
+ Add `commonVariantSitesIndex` input to CollectAllelicCounts.
+ Centrifuge: Fix issue where Centrifuge could not locate index files.
+ Increase default memory of BWA mem to 32G (was 16G).
+ Add `memory` input to fastqc task.
+ Centrifuge: Fix issue where centrifuge would fail on incorrect paths.
+ Added GATK CNV calling tasks:
    + AnnotateIntervals
    + CallCopyRatioSegments
    + CollectAllelicCounts
    + CollectReadCounts
    + CreateReadCountPanelOfNormals
    + DenoiseReadCounts
    + ModelSegments
    + PlotDenoisedCopyRatios
    + PlotModeledSegments
    + PreprocessIntervals
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
+ Picard's BedToIntervalList outputPath input is now
  optional (with a default of "regions.interval_list").
+ TALON: Fix SQLite error concerning database/disk space being full.
+ Update htseq to default image version 0.11.2.
+ Update biowdl-input-converter in common.wdl to version 0.2.1.
+ Update TALON section to now include the new annotation file output, and
  add config file creation to the TALON task.
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
+ Updated parameter_meta sections for Minimap2 and TranscriptClean to
  wdl-aid format.
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
+ Memory runtime attributes are now Strings indicating total memory, as
  opposed to Ints indicating memory per core.
+ Memory inputs for most tasks are now Strings, remaining Int memory inputs
  are renamed to "memoryGb".
+ Use the biowdl-input-converter container for JsonToYaml, to reduce the
  amount of containers needed.
+ Add biowdl-input-converter and remove SampleConfigToSampleReadgroupLists
  which it replaces.
+ GATK.GenotypeGVCFs: Increased memoryMultiplier from 2.0 to 3.0.
+ Minimap2: Add -k option to minimap2 mapping.
+ Added bwakit task.
+ Minimap2: Add the option for --MD tag.
+ TALON: Update average memory needs for main TALON process.

version 1.0.0
---------------------------
+ Common: Add "SampleConfigToSampleReadgroupLists" task.
+ MultiQC: the "interactive" input is now set to true by default.
+ Removed deprecated tasks: bioconda.installPrefix, mergecounts.MergeCounts
+ GATK.BaseRecalibrator: "knownIndelsSitesVCFs"
  and "knownIndelsSitesVCFIndexes" are no longer optional, but
  now have a default of "[]".
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
+ GATK: Add CombineVariants task that allows, e.g., to merge VCFs
  from different callers.
+ Mutect2: Add GATK tasks related to variant
  filtering (LearnReadOrientationModel, MergeStats, GetPileupSummaries,
  CalculateContamination and FilterMutectCalls).
+ Mutect2: Add "--germline-resource" and "--f1r2-tar-gz" inputs, requiring
  an update to GATK 4.1.2.0.
+ Mutect2: Add necessary missing index attribute for panel of normals.
+ MultiQC: Add memory variable to multiqc task.
+ GATK: SplitNCigarReads, BaseRecalibration and ApplyBQSR do no longer need
  regions files as required inputs.
+ VarDict: Add user definable flags (-M, -A, -Q, -d, -v, -f) to the paired
  VCF filtering script.
+ Cutadapt: If the output is a gzipped file, compress with
  level 1 (instead of default 6).
+ Cutadapt: Fix issues with read2output when using single-end reads.
+ Add feature type, idattr and additional attributes to htseq-count.
+ Added allow-contain option to bowtie.
+ Added a changelog to keep track of changes.
+ Added sortByName task in samtools to support more memory efficient
  execution of HTSeqCount.
+ Removed the bam index from HTSeqCount's inputs.
