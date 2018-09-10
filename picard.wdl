version 1.0

import "common.wdl"

task BedToIntervalList {
    input {
        String? preCommand
        File? picardJar

        File bedFile
        File dict
        String outputPath

        Int memory = 4
        Float memoryMultiplier = 3.0
    }

    String toolCommand = if defined(picardJar)
        then "java -Xmx" + memory + "G -jar " + picardJar
        else "picard -Xmx" + memory + "G"

    command {
        set -e -o pipefail
        mkdir -p $(dirname "~{outputPath}")
        ~{preCommand}
        ~{toolCommand} \
        BedToIntervalList \
        I=~{bedFile} \
        O=~{outputPath} \
        SD=~{dict}
    }

    output {
        File intervalList = outputPath
    }

    runtime {
        memory: ceil(memory * memoryMultiplier)
    }
}

task CollectMultipleMetrics {
    input {
        String? preCommand
        IndexedBamFile bam
        Reference reference
        String basename

        Boolean collectAlignmentSummaryMetrics = true
        Boolean collectInsertSizeMetrics = true
        Boolean qualityScoreDistribution = true
        Boolean meanQualityByCycle = true
        Boolean collectBaseDistributionByCycle = true
        Boolean collectGcBiasMetrics = true
        #FIXME: Boolean rnaSeqMetrics = false # There is a bug in picard https://github.com/broadinstitute/picard/issues/999
        Boolean collectSequencingArtifactMetrics = true
        Boolean collectQualityYieldMetrics = true

        String? picardJar

        Int memory = 4
        Float memoryMultiplier = 3.0
    }

    String toolCommand = if defined(picardJar)
        then "java -Xmx" + memory + "G -jar " + picardJar
        else "picard -Xmx" + memory + "G"

    command {
        set -e -o pipefail
        mkdir -p $(dirname "~{basename}")
        ~{preCommand}
        ~{toolCommand} \
        CollectMultipleMetrics \
        I=~{bam.file} \
        R=~{reference.fasta} \
        O=~{basename} \
        PROGRAM=null \
        ~{true="PROGRAM=CollectAlignmentSummaryMetrics" false="" collectAlignmentSummaryMetrics} \
        ~{true="PROGRAM=CollectInsertSizeMetrics" false="" collectInsertSizeMetrics} \
        ~{true="PROGRAM=QualityScoreDistribution" false="" qualityScoreDistribution} \
        ~{true="PROGRAM=MeanQualityByCycle" false="" meanQualityByCycle} \
        ~{true="PROGRAM=CollectBaseDistributionByCycle" false="" collectBaseDistributionByCycle} \
        ~{true="PROGRAM=CollectGcBiasMetrics" false="" collectGcBiasMetrics} \
        ~{true="PROGRAM=CollectSequencingArtifactMetrics" false=""
            collectSequencingArtifactMetrics} \
        ~{true="PROGRAM=CollectQualityYieldMetrics" false="" collectQualityYieldMetrics}
    }

    output {
        File alignmentSummary = basename + ".alignment_summary_metrics"
        File baitBiasDetail = basename + ".bait_bias_detail_metrics"
        File baitBiasSummary = basename + ".bait_bias_summary_metrics"
        File baseDistributionByCycle = basename + ".base_distribution_by_cycle_metrics"
        File baseDistributionByCyclePdf = basename + ".base_distribution_by_cycle.pdf"
        File errorSummary = basename + ".error_summary_metrics"
        File gcBiasDetail = basename + ".gc_bias.detail_metrics"
        File gcBiasPdf = basename + ".gc_bias.pdf"
        File gcBiasSummary = basename + ".gc_bias.summary_metrics"
        File? insertSizeHistogramPdf = basename + ".insert_size_histogram.pdf"
        File? insertSize = basename + ".insert_size_metrics"
        File preAdapterDetail = basename + ".pre_adapter_detail_metrics"
        File preAdapterSummary = basename + ".pre_adapter_summary_metrics"
        File qualityByCycle = basename + ".quality_by_cycle_metrics"
        File qualityByCyclePdf = basename + ".quality_by_cycle.pdf"
        File qualityDistribution = basename + ".quality_distribution_metrics"
        File qualityDistributionPdf = basename + ".quality_distribution.pdf"
        File qualityYield = basename + ".quality_yield_metrics"
    }

    runtime {
        memory: ceil(memory * memoryMultiplier)
    }
}

task CollectRnaSeqMetrics {
    input {
        String? preCommand
        IndexedBamFile bam
        File refRefflat
        String basename
        String strandSpecificity = "NONE"

        String? picardJar

        Int memory = 4
        Float memoryMultiplier = 3.0
    }

    String toolCommand = if defined(picardJar)
        then "java -Xmx" + memory + "G -jar " + picardJar
        else "picard -Xmx" + memory + "G"

    command {
        set -e -o pipefail
        mkdir -p $(dirname "~{basename}")
        ~{preCommand}
        ~{toolCommand} \
        CollectRnaSeqMetrics \
        I=~{bam.file} \
        O=~{basename}.RNA_Metrics \
        CHART_OUTPUT=~{basename}.RNA_Metrics.pdf \
        STRAND_SPECIFICITY=~{strandSpecificity} \
        REF_FLAT=~{refRefflat}
    }

    output {
        File? chart = basename + ".RNA_Metrics.pdf"
        File metrics = basename + ".RNA_Metrics"
    }

    runtime {
        memory: ceil(memory * memoryMultiplier)
    }
}

task CollectTargetedPcrMetrics {
    input {
        String? preCommand
        IndexedBamFile bam
        Reference reference
        File ampliconIntervals
        Array[File]+ targetIntervals
        String basename

        String? picardJar

        Int memory = 4
        Float memoryMultiplier = 3.0
    }

    String toolCommand = if defined(picardJar)
        then "java -Xmx" + memory + "G -jar " + picardJar
        else "picard -Xmx" + memory + "G"

    command {
        set -e -o pipefail
        mkdir -p $(dirname "~{basename}")
        ~{preCommand}
        ~{toolCommand} \
        CollectTargetedPcrMetrics \
        I=~{bam.file} \
        R=~{reference.fasta} \
        AMPLICON_INTERVALS=~{ampliconIntervals} \
        TARGET_INTERVALS=~{sep=" TARGET_INTERVALS=" targetIntervals} \
        O=~{basename}.targetPcrMetrics \
        PER_BASE_COVERAGE=~{basename}.targetPcrPerBaseCoverage \
        PER_TARGET_COVERAGE=~{basename}.targetPcrPerTargetCoverage
    }

    output {
        File perTargetCoverage = basename + ".targetPcrPerTargetCoverage"
        File perBaseCoverage = basename + ".targetPcrPerBaseCoverage"
        File metrics = basename + ".targetPcrMetrics"
    }

    runtime {
        memory: ceil(memory * memoryMultiplier)
    }
}

# Combine multiple recalibrated BAM files from scattered ApplyRecalibration runs
task GatherBamFiles {
    input {
        String? preCommand
        Array[IndexedBamFile]+ inputBams
        String outputBamPath
        String? picardJar

        Int memory = 4
        Float memoryMultiplier = 3.0
    }

    String toolCommand = if defined(picardJar)
        then "java -Xmx" + memory + "G -jar " + picardJar
        else "picard -Xmx" + memory + "G"

    command {
        set -e -o pipefail
        ~{preCommand}
        ~{toolCommand} \
        GatherBamFiles \
        INPUT=~{sep=' INPUT=' inputBams.file} \
        OUTPUT=~{outputBamPath} \
        CREATE_INDEX=true \
        CREATE_MD5_FILE=true
    }

    output {
        IndexedBamFile outputBam = object {
          file: outputBamPath,
          index: sub(outputBamPath, ".bam$", ".bai"),
          md5: outputBamPath + ".md5"
        }
    }

    runtime {
        memory: ceil(memory * memoryMultiplier)
    }
}

# Mark duplicate reads to avoid counting non-independent observations
task MarkDuplicates {
    input {
        String? preCommand
        Array[IndexedBamFile] inputBams
        String outputBamPath
        String metricsPath
        String? picardJar

        Int memory = 4
        Float memoryMultiplier = 3.0

        # The program default for READ_NAME_REGEX is appropriate in nearly every case.
        # Sometimes we wish to supply "null" in order to turn off optical duplicate detection
        # This can be desirable if you don't mind the estimated library size being wrong and optical duplicate detection is taking >7 days and failing
        String? read_name_regex
    }

    # Task is assuming query-sorted input so that the Secondary and Supplementary reads get marked correctly
    # This works because the output of BWA is query-grouped and therefore, so is the output of MergeBamAlignment.
    # While query-grouped isn't actually query-sorted, it's good enough for MarkDuplicates with ASSUME_SORT_ORDER="queryname"

    String toolCommand = if defined(picardJar)
        then "java -Xmx" + memory + "G -jar " + picardJar
        else "picard -Xmx" + memory + "G"

    command {
        set -e -o pipefail
        ~{preCommand}
        mkdir -p $(dirname ~{outputBamPath})
        ~{toolCommand} \
        MarkDuplicates \
        INPUT=~{sep=' INPUT=' inputBams.file} \
        OUTPUT=~{outputBamPath} \
        METRICS_FILE=~{metricsPath} \
        VALIDATION_STRINGENCY=SILENT \
        ~{"READ_NAME_REGEX=" + read_name_regex} \
        OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
        CLEAR_DT="false" \
        CREATE_INDEX=true \
        ADD_PG_TAG_TO_READS=false \
        CREATE_MD5_FILE=true
    }

    output {
        IndexedBamFile outputBam = object {
          file: outputBamPath,
          index: sub(outputBamPath, ".bam$", ".bai"),
          md5: outputBamPath + ".md5"
        }
        File metricsFile = metricsPath
    }

    runtime {
        memory: ceil(memory * memoryMultiplier)
    }
}

# Combine multiple VCFs or GVCFs from scattered HaplotypeCaller runs
task MergeVCFs {
    input {
        String? preCommand
        Array[File] inputVCFs
        Array[File] inputVCFsIndexes
        String outputVcfPath
        Int? compressionLevel
        String? picardJar

        Int memory = 4
        Float memoryMultiplier = 3.0
    }

    # Using MergeVcfs instead of GatherVcfs so we can create indices
    # See https://github.com/broadinstitute/picard/issues/789 for relevant GatherVcfs ticket

    String toolCommand = if defined(picardJar)
        then "java -Xmx" + memory + "G -jar " + picardJar
        else "picard -Xmx" + memory + "G"

    command {
        set -e -o pipefail
        ~{preCommand}
        ~{toolCommand} \
        MergeVcfs \
        INPUT=~{sep=' INPUT=' inputVCFs} \
        OUTPUT=~{outputVcfPath}
    }

    output {
        IndexedVcfFile outputVcf = object {
          file: outputVcfPath,
          index: outputVcfPath + ".tbi"
        }
    }

    runtime {
        memory: ceil(memory * memoryMultiplier)
    }
}

task SamToFastq {
    input {
        String? preCommand
        IndexedBamFile inputBam
        String outputRead1
        String? outputRead2
        String? outputUnpaired

        String? picardJar
        Int memory = 16 # High memory default to avoid crashes.
        Float memoryMultiplier = 3.0
    }

    String toolCommand = if defined(picardJar)
    then "java -Xmx" + memory + "G -jar " + picardJar
    else "picard -Xmx" + memory + "G"

    command {
        set -e -o pipefail
        ~{preCommand}
        ~{toolCommand} \
        SamToFastq \
        I=~{inputBam.file} \
        ~{"FASTQ=" + outputRead1} \
        ~{"SECOND_END_FASTQ=" + outputRead2} \
        ~{"UNPAIRED_FASTQ=" + outputUnpaired}
    }

    output {
        File read1 = outputRead1
        File? read2 = outputRead2
        File? unpairedRead = outputUnpaired
    }

    runtime {
        memory: ceil(memory * memoryMultiplier)
    }
}

task ScatterIntervalList {
    input {
        String? preCommand
        File interval_list
        Int scatter_count
        String? picardJar

        Int memory = 4
        Float memoryMultiplier = 3.0
    }

    String toolCommand = if defined(picardJar)
        then "java -Xmx" + memory + "G -jar " + picardJar
        else "picard -Xmx" + memory + "G"

    command {
        set -e -o pipefail
        ~{preCommand}
        mkdir scatter_list
        ~{toolCommand} \
        IntervalListTools \
        SCATTER_COUNT=~{scatter_count} \
        SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
        UNIQUE=true \
        SORT=true \
        INPUT=~{interval_list} \
        OUTPUT=scatter_list
    }

    output {
        Array[File] out = glob("scatter_list/*/*.interval_list")
        Int interval_count = read_int(stdout())
    }

    runtime {
        memory: ceil(memory * memoryMultiplier)
    }
}

task SortVcf {
    input {
        String? preCommand
        String? picardJar

        Array[File]+ vcfFiles
        String outputVcfPath
        File? dict

        Int memory = 4
        Float memoryMultiplier = 3.0
        }

        String toolCommand = if defined(picardJar)
            then "java -Xmx" + memory + "G -jar " + picardJar
            else "picard -Xmx" + memory + "G"

    command {
        set -e -o pipefail
        ~{preCommand}
        ~{toolCommand} \
        SortVcf \
        I=~{sep=" I=" vcfFiles} \
        ~{"SEQUENCE_DICTIONARY=" + dict} \
        O=~{outputVcfPath}
    }

    output {
        IndexedVcfFile outputVcf = object {
          file: outputVcfPath,
          index: outputVcfPath + ".tbi"
        }
    }

    runtime {
        memory: ceil(memory * memoryMultiplier)
    }
}