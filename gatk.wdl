version 1.0

# Copyright (c) 2018 Leiden University Medical Center
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

task AnnotateIntervals {
    input {
        File referenceFasta
        File referenceFastaDict
        File referenceFastaFai
        String annotatedIntervalsPath = "intervals.annotated.tsv"
        File intervals
        String intervalMergingRule = "OVERLAPPING_ONLY"
        Int featureQueryLookahead = 1000000

        File? mappabilityTrack
        File? segmentalDuplicationTrack

        String javaXmx = "2G"
        String memory = "3GiB"
        Int timeMinutes = 5
        String dockerImage = "quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{annotatedIntervalsPath})"
        gatk --java-options '-Xmx~{javaXmx} -XX:ParallelGCThreads=1' \
        AnnotateIntervals \
        -R ~{referenceFasta} \
        -L ~{intervals} \
        ~{"--mappability-track  " + mappabilityTrack} \
        ~{"--segmental-duplication-track " + segmentalDuplicationTrack} \
        --feature-query-lookahead ~{featureQueryLookahead} \
        --interval-merging-rule ~{intervalMergingRule} \
        -O ~{annotatedIntervalsPath}
    }

    output {
        File annotatedIntervals = annotatedIntervalsPath
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        referenceFasta: {description: "The reference fasta file.", category: "required"}
        referenceFastaDict: {description: "The sequence dictionary associated with the reference fasta file.", category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        annotatedIntervalsPath: {description: "The location the output should be written to.", category: "advanced"}
        intervals: {description: "An interval list describinig the intervals to annotate.", category: "required"}
        intervalMergingRule: {description: "Equivalent to gatk AnnotateIntervals' `--interval-merging-rule` option.", category: "advanced"}
        mappabilityTrack: {description: "Equivalent to gatk AnnotateIntervals' `--mappability-track` option.", category: "common"}
        segmentalDuplicationTrack: {description: "Equivalent to gatk AnnotateIntervals' `--segmenta-duplicarion-track` option.", category: "common"}
        featureQueryLookahead: {description: "Equivalent to gatk AnnotateIntervals' `--feature-query-lookahead` option.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        annotatedIntervals: {description: "This is a tab-separated values (TSV) file with a SAM-style header containing a sequence dictionary, a row specifying the column headers for the contained annotations, and the corresponding entry rows."}
    }
}

# Apply Base Quality Score Recalibration (BQSR) model.
task ApplyBQSR {
    input {
        File inputBam
        File inputBamIndex
        String outputBamPath
        File recalibrationReport
        Array[File] sequenceGroupInterval = []
        File referenceFasta
        File referenceFastaDict
        File referenceFastaFai

        Int javaXmxMb = 2048
        Int memoryMb = javaXmxMb + 512
        # This will likely be used with intervals, as such size based
        # estimation can't be used.
        Int timeMinutes = 120
        String dockerImage = "quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputBamPath})"
        gatk --java-options '-Xmx~{javaXmxMb}M -XX:ParallelGCThreads=1' \
        ApplyBQSR \
        --create-output-bam-md5 \
        --add-output-sam-program-record \
        -R ~{referenceFasta} \
        -I ~{inputBam} \
        --use-original-qualities \
        -O ~{outputBamPath} \
        -bqsr ~{recalibrationReport} \
        --static-quantized-quals 10 \
        --static-quantized-quals 20 \
        --static-quantized-quals 30 \
        ~{true="-L" false="" length(sequenceGroupInterval) > 0} ~{sep=' -L ' sequenceGroupInterval}
    }

    output {
        File recalibratedBam = outputBamPath
        File recalibratedBamIndex = sub(outputBamPath, "\.bam$", ".bai")
        File recalibratedBamMd5 = outputBamPath + ".md5"
    }

    runtime {
        memory: "~{memoryMb}MiB"
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        inputBam: {description: "The BAM file which should be recalibrated.", category: "required"}
        inputBamIndex: {description: "The input BAM file's index.", category: "required"}
        outputBamPath: {description: "The location the resulting BAM file should be written.", category: "required"}
        recalibrationReport: {description: "The BQSR report the be used for recalibration.", category: "required"}
        sequenceGroupInterval: {description: "Bed files describing the regions to operate on.", category: "advanced"}
        referenceFasta: {description: "The reference fasta file which was also used for mapping.", category: "required"}
        referenceFastaDict: {description: "The sequence dictionary associated with the reference fasta file.", category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        javaXmxMb: {description: "The maximum memory available to the program in megabytes. Should be lower than `memoryMb` to accommodate JVM overhead.", category: "advanced"}
        memoryMb: {description: "The amount of memory this job will use in megabytes.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        recalibratedBam: {description: "A BAM file containing the recalibrated read data."}
        recalibratedBamIndex: {description: "Index of recalibrated BAM file."}
        recalibratedBamMd5: {description: "MD5 of recalibrated BAM file."}
    }
}

# Generate Base Quality Score Recalibration (BQSR) model.
task BaseRecalibrator {
    input {
        File inputBam
        File inputBamIndex
        String recalibrationReportPath
        Array[File] sequenceGroupInterval = []
        Array[File] knownIndelsSitesVCFs = []
        Array[File] knownIndelsSitesVCFIndexes = []
        File referenceFasta
        File referenceFastaDict
        File referenceFastaFai

        File? dbsnpVCF
        File? dbsnpVCFIndex

        Int javaXmxMb = 1024
        Int memoryMb = javaXmxMb + 512
        Int timeMinutes = 120 # This will likely be used with intervals, as such size based estimation can't be used.
        String dockerImage = "quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{recalibrationReportPath})"
        gatk --java-options '-Xmx~{javaXmxMb}M -XX:ParallelGCThreads=1' \
        BaseRecalibrator \
        -R ~{referenceFasta} \
        -I ~{inputBam} \
        --use-original-qualities \
        -O ~{recalibrationReportPath} \
        ~{true="--known-sites" false="" length(knownIndelsSitesVCFs) > 0} ~{sep=" --known-sites " knownIndelsSitesVCFs} \
        ~{"--known-sites " + dbsnpVCF} \
        ~{true="-L" false="" length(sequenceGroupInterval) > 0} ~{sep=' -L ' sequenceGroupInterval}
    }

    output {
        File recalibrationReport = recalibrationReportPath
    }

    runtime {
        memory: "~{memoryMb}MiB"
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        inputBam: {description: "The BAM file to generate a BQSR report for.", category: "required"}
        inputBamIndex: {description: "The index of the input BAM file.", category: "required"}
        recalibrationReportPath: {description: "The location to write the BQSR report to.", category: "required"}
        sequenceGroupInterval: {description: "Bed files describing the regions to operate on.", category: "advanced"}
        knownIndelsSitesVCFs: {description: "VCF files with known indels.", category: "advanced"}
        knownIndelsSitesVCFIndexes: {description: "The indexed for the known variant VCFs.", category: "advanced"}
        referenceFasta: {description: "The reference fasta file which was also used for mapping.", category: "required"}
        referenceFastaDict: {description: "The sequence dictionary associated with the reference fasta file.", category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        dbsnpVCF: {description: "A dbSNP VCF.", category: "common"}
        dbsnpVCFIndex: {description: "The index for the dbSNP VCF.", category: "common"}
        javaXmxMb: {description: "The maximum memory available to the program in megabytes. Should be lower than `memoryMb` to accommodate JVM overhead.", category: "advanced"}
        memoryMb: {description: "The amount of memory this job will use in megabytes.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        recalibrationReport: {description: "A GATK Report file with many tables."}
    }
}

task CalculateContamination {
    input {
        File tumorPileups

        File? normalPileups

        String javaXmx = "12G"
        String memory = "13GiB"
        Int timeMinutes = 180
        String dockerImage = "quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0"
    }

    command {
        set -e
        gatk --java-options '-Xmx~{javaXmx} -XX:ParallelGCThreads=1' \
        CalculateContamination \
        -I ~{tumorPileups} \
        ~{"-matched " + normalPileups} \
        -O "contamination.table" \
        --tumor-segmentation "segments.table"
    }

    output {
        File contaminationTable = "contamination.table"
        File mafTumorSegments = "segments.table"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        tumorPileups: {description: "The pileup summary of a tumor/case sample.", category: "required"}
        normalPileups: {description: "The pileup summary of the normal/control sample.", category: "common"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        contaminationTable: {description: "Table with fractions of reads from cross-sample contamination."}
        mafTumorSegments: {description: "Tumor segments table."}
    }
}

task CallCopyRatioSegments {
    input {
        String outputPrefix
        File copyRatioSegments

        String javaXmx = "2G"
        String memory = "3GiB"
        Int timeMinutes = 2
        String dockerImage = "quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPrefix})"
        gatk --java-options '-Xmx~{javaXmx} -XX:ParallelGCThreads=1' \
        CallCopyRatioSegments \
        -I ~{copyRatioSegments} \
        -O ~{outputPrefix}.called.seg
    }

    output {
        File calledSegments = outputPrefix + ".called.seg"
        File calledSegmentsIgv = outputPrefix + ".called.igv.seg"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        outputPrefix: {description: "The prefix for the output files.", category: "required"}
        copyRatioSegments: {description: "The copy ratios file generated by gatk ModelSegments.", category: "required"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        calledSegments: {description: "This is a tab-separated values (TSV) file with a SAM-style header containing a read group sample name, a sequence dictionary, a row specifying the column headers contained in CalledCopyRatioSegmentCollection.CalledCopyRatioSegmentTableColumn, and the corresponding entry rows."}
        calledSegmentsIgv: {description: "This is a tab-separated values (TSV) file with CBS-format column headers and the corresponding entry rows that can be plotted using IGV."}
    }
}

task CollectAllelicCounts {
    input {
        String allelicCountsPath = "allelic_counts.tsv"
        File commonVariantSites
        File inputBam
        File inputBamIndex
        File referenceFasta
        File referenceFastaDict
        File referenceFastaFai

        File? commonVariantSitesIndex

        String javaXmx = "10G"
        String memory = "11GiB"
        Int timeMinutes = 120
        String dockerImage = "quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{allelicCountsPath})"
        gatk --java-options '-Xmx~{javaXmx} -XX:ParallelGCThreads=1' \
        CollectAllelicCounts \
        -L ~{commonVariantSites} \
        -I ~{inputBam} \
        -R ~{referenceFasta} \
        -O ~{allelicCountsPath}
    }

    output {
        File allelicCounts = allelicCountsPath
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        allelicCountsPath: {description: "The path the output should be written to.", category: "advanced"}
        commonVariantSites: {description: "Interval list or vcf of common variant sites (to retrieve the allelic counts for).", category: "required"}
        inputBam: {description: "The BAM file to generate counts for.", category: "required"}
        inputBamIndex: {description: "The index of the input BAM file.", category: "required"}
        referenceFasta: {description: "The reference fasta file.", category: "required"}
        referenceFastaDict: {description: "The sequence dictionary associated with the reference fasta file.", category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        commonVariantSitesIndex: {description: "The index for commonVariantSites.", category: "common"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        allelicCounts: {description: "This is a tab-separated values (TSV) file with a SAM-style header containing a read group sample name, a sequence dictionary, a row specifying the column headers contained in AllelicCountCollection.AllelicCountTableColumn, and the corresponding entry rows."}
    }
}

task CollectReadCounts {
    input {
        String countsPath = "readcounts.hdf5"
        File intervals
        File inputBam
        File inputBamIndex
        File referenceFasta
        File referenceFastaDict
        File referenceFastaFai
        String intervalMergingRule = "OVERLAPPING_ONLY"

        String javaXmx = "7G"
        String memory = "8GiB"
        Int timeMinutes = 1 + ceil(size(inputBam, "G") * 5)
        String dockerImage = "quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{countsPath})"
        gatk --java-options '-Xmx~{javaXmx} -XX:ParallelGCThreads=1' \
        CollectReadCounts \
        -L ~{intervals} \
        -I ~{inputBam} \
        -R ~{referenceFasta} \
        --format HDF5 \
        --interval-merging-rule ~{intervalMergingRule} \
        -O ~{countsPath}
    }

    output {
        File counts = countsPath
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        countsPath: {description: "The location the output should be written to.", category: "advanced"}
        intervals: {description: "The intervals to collect counts for.", category: "required"}
        inputBam: {description: "The BAM file to determine the coverage for.", category: "required"}
        inputBamIndex: {description: "The input BAM file's index.", category: "required"}
        referenceFasta: {description: "The reference fasta file.", category: "required"}
        referenceFastaDict: {description: "The sequence dictionary associated with the reference fasta file.", category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        intervalMergingRule: {description: "Equivalent to gatk CollectReadCounts' `--interval-merging-rule` option.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        counts: {description: "Read counts at specified intervals."}
    }
}

task CombineGVCFs {
    input {
        Array[File]+ gvcfFiles
        Array[File]+ gvcfFilesIndex
        Array[File] intervals = []
        String outputPath
        File referenceFasta
        File referenceFastaDict
        File referenceFastaFai

        String javaXmx = "4G"
        String memory = "5GiB"
        Int timeMinutes = 1 + ceil(size(gvcfFiles, "G") * 8)
        String dockerImage = "quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPath})"
        gatk --java-options '-Xmx~{javaXmx} -XX:ParallelGCThreads=1' \
        CombineGVCFs \
        -R ~{referenceFasta} \
        -O ~{outputPath} \
        -V ~{sep=' -V ' gvcfFiles} \
        ~{true='-L' false='' length(intervals) > 0} ~{sep=' -L ' intervals}
    }

    output {
        File outputVcf = outputPath
        File outputVcfIndex = outputPath + ".tbi"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        gvcfFiles: {description: "The GVCF files to be combined.", category: "required"}
        gvcfFilesIndex: {description: "The indexes for the GVCF files.", caregory: "required"}
        intervals: {description: "Bed files or interval lists describing the regions to operate on.", category: "advanced"}
        outputPath: {description: "The location the combined GVCF should be written to.", category: "required"}
        referenceFasta: {description: "The reference fasta file which was also used for mapping.", category: "required"}
        referenceFastaDict: {description: "The sequence dictionary associated with the reference fasta file.", category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        outputVcf: {description: "A combined multi-sample gVCF."}
        outputVcfIndex: {description: "Index of the output file."}
    }
}

task CombineVariants {
    input {
        File referenceFasta
        File referenceFastaFai
        File referenceFastaDict
        String genotypeMergeOption = "UNIQUIFY"
        String filteredRecordsMergeType = "KEEP_IF_ANY_UNFILTERED"
        Array[String]+ identifiers
        Array[File]+ variantVcfs # Follow "identifiers" array order.
        Array[File]+ variantIndexes
        String outputPath

        String javaXmx = "12G"
        String memory = "13GiB"
        Int timeMinutes = 180
        String dockerImage = "broadinstitute/gatk3:3.8-1"
    }

    command <<<
        set -e
        mkdir -p "$(dirname ~{outputPath})"
        # Build "-V:<ID> <file.vcf>" arguments according to IDs
        # and VCFs to merge.
        # Make sure commands are run in bash.
        V_args=$(bash -c '
        set -eu
        ids=(~{sep=" " identifiers})
        vars=(~{sep=" " variantVcfs})
        for (( i = 0; i < ${#ids[@]}; ++i ))
        do
            printf -- "-V:%s %s " "${ids[i]}" "${vars[i]}"
        done
        ')
        java -Xmx~{javaXmx} -XX:ParallelGCThreads=1 -jar /usr/GenomeAnalysisTK.jar \
        -T CombineVariants \
        -R ~{referenceFasta} \
        --genotypemergeoption ~{genotypeMergeOption} \
        --filteredrecordsmergetype ~{filteredRecordsMergeType} \
        --out ~{outputPath} \
        $V_args
    >>>

    output {
        File combinedVcf = outputPath
        File combinedVcfIndex = outputPath + ".tbi"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        referenceFasta: {description: "The reference fasta file which was also used for mapping.", category: "required"}
        referenceFastaDict: {description: "The sequence dictionary associated with the reference fasta file.", category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        genotypeMergeOption: {description: "Equivalent to CombineVariants' `--genotypemergeoption` option.", category: "advanced"}
        filteredRecordsMergeType: {description: "Equivalent to CombineVariants' `--filteredrecordsmergetype` option.", category: "advanced"}
        identifiers: {description: "The sample identifiers in the same order as variantVcfs.", category: "required"}
        variantVcfs: {description: "The input VCF files in the same order as identifiers.", category: "required"}
        variantIndexes: {description: "The indexes of the input VCF files.", category: "required"}
        outputPath: {description: "The location the output should be written to.", category: "required"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        combinedVcf: {description: "Combined VCF file."}
        combinedVcfIndex: {description: "Index of combined VCF file."}
    }
}

task CreateReadCountPanelOfNormals {
    input {
        String PONpath = "PON.hdf5"
        Array[File]+ readCountsFiles

        File? annotatedIntervals

        String javaXmx = "7G"
        String memory = "8GiB"
        Int timeMinutes = 5
        # The biocontainer causes a spark related error for some reason.
        String dockerImage = "broadinstitute/gatk:4.1.8.0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{PONpath})"
        gatk --java-options '-Xmx~{javaXmx} -XX:ParallelGCThreads=1' \
        CreateReadCountPanelOfNormals \
        -I ~{sep=" -I " readCountsFiles} \
        ~{"--annotated-intervals " + annotatedIntervals} \
        -O ~{PONpath}
    }

    output {
        File PON = PONpath
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        PONpath: {description: "The location the PON should be written to.", category: "common"}
        readCountsFiles: {description: "The read counts files as generated by CollectReadCounts.", category: "required"}
        annotatedIntervals: {description: "An annotation set of intervals as generated by AnnotateIntervals. If provided, explicit GC correction will be performed.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        PON: {description: "Panel-of-normals file."}
    }
}

task DenoiseReadCounts {
    input {
        File readCounts
        String outputPrefix

        File? PON
        File? annotatedIntervals

        String javaXmx = "4G"
        String memory = "5GiB"
        Int timeMinutes = 5
        String dockerImage = "quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPrefix})"
        gatk --java-options '-Xmx~{javaXmx} -XX:ParallelGCThreads=1' \
        DenoiseReadCounts \
        -I ~{readCounts} \
        ~{"--count-panel-of-normals " + PON} \
        ~{"--annotated-intervals " + annotatedIntervals} \
        --standardized-copy-ratios ~{outputPrefix}.standardizedCR.tsv \
        --denoised-copy-ratios ~{outputPrefix}.denoisedCR.tsv
    }

    output {
        File standardizedCopyRatios = outputPrefix + ".standardizedCR.tsv"
        File denoisedCopyRatios = outputPrefix + ".denoisedCR.tsv"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        readCounts: {description: "The read counts file as generated by CollectReadCounts.", category: "required"}
        outputPrefix: {description: "The prefix for the output files.", category: "required"}
        PON: {description: "A panel of normals as generated by CreateReadCountPanelOfNormals.", category: "advanced"}
        annotatedIntervals: {description: "An annotated set of intervals as generated by AnnotateIntervals. Will be ignored if PON is provided.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        standardizedCopyRatios: {description: "This is a tab-separated values (TSV) file with a SAM-style header containing a read group sample name, a sequence dictionary, a row specifying the column headers contained in CopyRatioCollection.CopyRatioTableColumn, and the corresponding entry rows."}
        denoisedCopyRatios: {description: "This is a tab-separated values (TSV) file with a SAM-style header containing a read group sample name, a sequence dictionary, a row specifying the column headers contained in CopyRatioCollection.CopyRatioTableColumn, and the corresponding entry rows."}
    }
}

task FilterMutectCalls {
    input {
        File referenceFasta
        File referenceFastaFai
        File referenceFastaDict
        File unfilteredVcf
        File unfilteredVcfIndex
        String outputVcf
        Int uniqueAltReadCount = 4
        File mutect2Stats

        File? contaminationTable
        File? mafTumorSegments
        File? artifactPriors

        String javaXmx = "12G"
        String memory = "13GiB"
        Int timeMinutes = 60
        String dockerImage = "quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputVcf})"
        gatk --java-options '-Xmx~{javaXmx} -XX:ParallelGCThreads=1' \
        FilterMutectCalls \
        -R ~{referenceFasta} \
        -V ~{unfilteredVcf} \
        -O ~{outputVcf} \
        ~{"--contamination-table " + contaminationTable} \
        ~{"--tumor-segmentation " + mafTumorSegments} \
        ~{"--ob-priors " + artifactPriors} \
        ~{"--unique-alt-read-count " + uniqueAltReadCount} \
        ~{"-stats " + mutect2Stats} \
        --filtering-stats "filtering.stats" \
        --showHidden
    }

    output {
        File filteredVcf = outputVcf
        File filteredVcfIndex = outputVcf + ".tbi"
        File filteringStats = "filtering.stats"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        referenceFasta: {description: "The reference fasta file which was also used for mapping.", category: "required"}
        referenceFastaDict: {description: "The sequence dictionary associated with the reference fasta file.", category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        unfilteredVcf: {description: "An unfiltered VCF file as produced by Mutect2.", category: "required"}
        unfilteredVcfIndex: {description: "The index of the unfiltered VCF file.", category: "required"}
        outputVcf: {description: "The location the filtered VCF file should be written.", category: "required"}
        uniqueAltReadCount: {description: "Equivalent to FilterMutectCalls' `--unique-alt-read-count` option.", category: "advanced"}
        mutect2Stats: {description: "Equivalent to FilterMutectCalls' `-stats` option.", category: "advanced"}
        contaminationTable: {description: "Equivalent to FilterMutectCalls' `--contamination-table` option.", category: "advanced"}
        mafTumorSegments: {description: "Equivalent to FilterMutectCalls' `--tumor-segmentation` option.", category: "advanced"}
        artifactPriors: {description: "Equivalent to FilterMutectCalls' `--ob-priors` option.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        filteredVcf: {description: "VCF file with filtered variants from a Mutect2 VCF callset."}
        filteredVcfIndex: {description: "Index of output VCF file."}
        filteringStats: {description: "The output filtering stats file."}
    }
}

# Combine multiple recalibration tables from scattered BaseRecalibrator runs.
task GatherBqsrReports {
    input {
        Array[File] inputBQSRreports
        String outputReportPath

        Int javaXmxMb = 256
        Int memoryMb = 256 + javaXmxMb
        Int timeMinutes = 1
        String dockerImage = "quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputReportPath})"
        gatk --java-options '-Xmx~{javaXmxMb}M -XX:ParallelGCThreads=1' \
        GatherBQSRReports \
        -I ~{sep=' -I ' inputBQSRreports} \
        -O ~{outputReportPath}
    }

    output {
        File outputBQSRreport = outputReportPath
    }

    runtime {
        memory: "~{memoryMb}MiB"
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        inputBQSRreports: {description: "The BQSR reports to be merged.", category: "required"}
        outputReportPath: {description: "The location of the combined BQSR report.", category: "required"}
        javaXmxMb: {description: "The maximum memory available to the program in megabytes. Should be lower than `memory` to accommodate JVM overhead.", category: "advanced"}
        memoryMb: {description: "The amount of memory this job will use in megabytes.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        outputBQSRreport: {description: "Single file with scattered BQSR recalibration reports gathered into one."}
    }
}

task GenomicsDBImport {
    input {
        Array[File] gvcfFiles
        Array[File] gvcfFilesIndex
        Array[File]+ intervals
        String genomicsDBWorkspacePath = "genomics_db"
        String genomicsDBTarFile = "genomics_db.tar.gz"

        String? tmpDir

        String javaXmx = "4G"
        String memory = "5GiB"
        Int timeMinutes = 180
        String dockerImage = "quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{genomicsDBWorkspacePath})"
        gatk --java-options '-Xmx~{javaXmx} -XX:ParallelGCThreads=1' \
        GenomicsDBImport \
        -V ~{sep=" -V " gvcfFiles} \
        --genomicsdb-workspace-path ~{genomicsDBWorkspacePath} \
        ~{"--tmp-dir " + tmpDir} \
        -L ~{sep=" -L " intervals}
        bash -c 'tar -cvzf ~{genomicsDBTarFile} ~{genomicsDBWorkspacePath}/*'
    }

    output {
        File genomicsDbTarArchive = genomicsDBTarFile
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        gvcfFiles: {description: "The gvcfFiles to be merged.", category: "required"}
        gvcfFilesIndex: {description: "Indexes for the gvcfFiles.", category: "required"}
        intervals: {description: "intervals over which to operate.", category: "required"}
        genomicsDBWorkspacePath: {description: "Where the genomicsDB files should be stored.", category: "advanced"}
        genomicsDBTarFile: {description: "Where the .tar file containing the genomicsDB should be stored.", category: "advanced"}
        tmpDir: {description: "Alternate temporary directory in case there is not enough space. Must be mounted when using containers.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        genomicsDbTarArchive: {description: "Imported VCFs to GenomicsDB file."}
    }
}

task GenotypeGVCFs {
    input {
        File gvcfFile
        File gvcfFileIndex
        String outputPath
        File referenceFasta
        File referenceFastaDict
        File referenceFastaFai
        Array[String] annotationGroups = ["StandardAnnotation"]

        Array[File]? intervals
        File? dbsnpVCF
        File? dbsnpVCFIndex
        File? pedigree

        String javaXmx = "6G"
        String memory = "7GiB"
        Int timeMinutes = 120 # This will likely be used with intervals, as such size based estimation can't be used.
        String dockerImage = "quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPath})"
        gatk --java-options '-Xmx~{javaXmx} -XX:ParallelGCThreads=1' \
        GenotypeGVCFs \
        -R ~{referenceFasta} \
        -O ~{outputPath} \
        ~{"-D " + dbsnpVCF} \
        ~{"--pedigree " + pedigree} \
        ~{true="-G" false="" length(annotationGroups) > 0} ~{sep=" -G " annotationGroups} \
        -V ~{gvcfFile} \
        ~{true="--only-output-calls-starting-in-intervals" false="" defined(intervals)} \
        ~{true="-L" false="" defined(intervals)} ~{sep=' -L ' intervals}
    }

    output {
        File outputVCF = outputPath
        File outputVCFIndex = outputPath + ".tbi"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        gvcfFile: {description: "The GVCF file to be genotyped.", category: "required"}
        gvcfFileIndex: {description: "The index of the input GVCF file.", category: "required"}
        outputPath: {description: "The location to write the output VCF file to.", category: "required"}
        referenceFasta: {description: "The reference fasta file which was also used for mapping.", category: "required"}
        referenceFastaDict: {description: "The sequence dictionary associated with the reference fasta file.", category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        annotationGroups: {description: "Which annotation groups will be used for the annotation.", category: "advanced"}
        intervals: {description: "Bed files or interval lists describing the regions to operate on.", category: "optional"}
        dbsnpVCF: {description: "A dbSNP VCF.", category: "common"}
        dbsnpVCFIndex: {description: "The index for the dbSNP VCF.", category: "common"}
        pedigree: {description: "Pedigree file for determining the population \"founders\".", category: "common"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        outputVCF: {description: "A final VCF in which all samples have been jointly genotyped. "}
        outputVCFIndex: {description: "Index of final VCF file."}
    }
}

task GetPileupSummaries {
    input {
        File sampleBam
        File sampleBamIndex
        File variantsForContamination
        File variantsForContaminationIndex
        File sitesForContamination
        File sitesForContaminationIndex
        String outputPrefix

        String javaXmx = "12G"
        String memory = "13GiB"
        Int timeMinutes = 120
        String dockerImage = "quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0"
    }

    command {
        set -e
        gatk --java-options '-Xmx~{javaXmx} -XX:ParallelGCThreads=1' \
        GetPileupSummaries \
        -I ~{sampleBam} \
        -V ~{variantsForContamination} \
        -L ~{sitesForContamination} \
        -O ~{outputPrefix + "-pileups.table"}
    }

    output {
        File pileups = outputPrefix + "-pileups.table"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        sampleBam: {description: "A BAM file for which a pileup should be created.", category: "required"}
        sampleBamIndex: {description: "The index of the input BAM file.", category: "required"}
        variantsForContamination: {description: "A VCF file with common variants.", category: "required"}
        variantsForContaminationIndex: {description: "The index for the common variants VCF file.", category: "required"}
        sitesForContamination: {description: "A bed file describing regions to operate on.", category: "required"}
        sitesForContaminationIndex: {description: "The index for the bed file.", category: "required"}
        outputPrefix: {description: "The prefix for the ouput.", category: "required"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        pileups: {description: "Pileup metrics for inferring contamination."}
    }
}


task HaplotypeCaller {
    input {
        Array[File]+ inputBams
        Array[File]+ inputBamsIndex
        String outputPath
        File referenceFasta
        File referenceFastaIndex
        File referenceFastaDict
        Boolean gvcf = false
        String emitRefConfidence = if gvcf then "GVCF" else "NONE"
        Boolean dontUseSoftClippedBases = false

        Array[File]+? intervalList
        Array[File]+? excludeIntervalList
        Float? contamination
        File? dbsnpVCF
        File? dbsnpVCFIndex
        File? pedigree
        Int? ploidy
        String? outputMode
        Float? standardMinConfidenceThresholdForCalling

        Int javaXmxMb = 4096
        # Memory increases with time used. 4G should cover most use cases.
        Int memoryMb = javaXmxMb + 512
        Int timeMinutes = 400 # This will likely be used with intervals, as such size based estimation can't be used.
        String dockerImage = "quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPath})"
        gatk --java-options '-Xmx~{javaXmxMb}M -XX:ParallelGCThreads=1' \
        HaplotypeCaller \
        -R ~{referenceFasta} \
        -O ~{outputPath} \
        -I ~{sep=" -I " inputBams} \
        ~{"--sample-ploidy " + ploidy} \
        ~{true="-L" false="" defined(intervalList)} ~{sep=' -L ' intervalList} \
        ~{true="-XL" false="" defined(excludeIntervalList)} ~{sep=' -XL ' excludeIntervalList} \
        ~{"-D " + dbsnpVCF} \
        ~{"--pedigree " + pedigree} \
        ~{"--contamination-fraction-per-sample-file " + contamination} \
        ~{"--output-mode " + outputMode} \
        --emit-ref-confidence ~{emitRefConfidence} \
        ~{true="--dont-use-soft-clipped-bases" false="" dontUseSoftClippedBases} \
        ~{"--standard-min-confidence-threshold-for-calling " + standardMinConfidenceThresholdForCalling}
    }

    output {
        File outputVCF = outputPath
        File outputVCFIndex = outputPath + ".tbi"
    }

    runtime {
        memory: "~{memoryMb}MiB"
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        inputBams: {description: "The BAM files on which to perform variant calling.", category: "required"}
        inputBamsIndex: {description: "The indexes for the input BAM files.", category: "required"}
        outputPath: {description: "The location to write the output to.", category: "required"}
        referenceFasta: {description: "The reference fasta file which was also used for mapping.", category: "required"}
        referenceFastaDict: {description: "The sequence dictionary associated with the reference fasta file.", category: "required"}
        referenceFastaIndex: {description: "The index for the reference fasta file.", category: "required"}
        gvcf: {description: "Whether the output should be a gvcf.", category: "common"}
        emitRefConfidence: {description: "Whether to include reference calls. Three modes: 'NONE', 'BP_RESOLUTION' and 'GVCF'.", category: "advanced"}
        dontUseSoftClippedBases: {description: "Do not use soft-clipped bases. Should be 'true' for RNA variant calling.", category: "common"}
        intervalList: {description: "Bed files or interval lists describing the regions to operate on.", category: "common"}
        excludeIntervalList: {description: "Bed files or interval lists describing the regions to NOT operate on.", category: "common"}
        contamination: {description: "Equivalent to HaplotypeCaller's `-contamination` option.", category: "advanced"}
        dbsnpVCF: {description: "A dbSNP VCF.", category: "common"}
        dbsnpVCFIndex: {description: "The index for the dbSNP VCF.", category: "common"}
        pedigree: {description: "Pedigree file for determining the population \"founders\".", category: "common"}
        ploidy: {description: "The ploidy with which the variants should be called.", category: "common"}
        outputMode: {description: "Specifies which type of calls we should output. Same as HaplotypeCaller's `--output-mode` option.", category: "advanced"}
        standardMinConfidenceThresholdForCalling: {description: "Confidence threshold used for calling variants.", category: "advanced"}
        javaXmxMb: {description: "The maximum memory available to the program in megabytes. Should be lower than `memoryMb` to accommodate JVM overhead.", category: "advanced"}
        memoryMb: {description: "The amount of memory this job will use in megabytes.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        outputVCF: {description: "Raw, unfiltered SNP and indel calls."}
        outputVCFIndex: {description: "Index of output VCF."}
    }
}

task LearnReadOrientationModel {
    input {
        Array[File]+ f1r2TarGz

        String javaXmx = "12G"
        String memory = "13GiB"
        Int timeMinutes = 120
        String dockerImage = "quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0"
    }

    command {
        set -e
        gatk --java-options '-Xmx~{javaXmx} -XX:ParallelGCThreads=1' \
        LearnReadOrientationModel \
        -I ~{sep=" -I " f1r2TarGz} \
        -O "artifact-priors.tar.gz"
    }

    output {
        File artifactPriorsTable = "artifact-priors.tar.gz"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        f1r2TarGz: {description: "A f1r2TarGz file outputed by mutect2.", category: "required"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        artifactPriorsTable: {description: "Maximum likelihood estimates of artifact prior probabilities in the orientation bias mixture model filter."}
    }
}

task MergeStats {
    input {
        Array[File]+ stats

        String javaXmx = "14G"
        String memory = "15GiB"
        Int timeMinutes = 30
        String dockerImage = "quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0"
    }

    command {
        set -e
        gatk --java-options '-Xmx~{javaXmx} -XX:ParallelGCThreads=1' \
        MergeMutectStats \
        -stats ~{sep=" -stats " stats} \
        -O "merged.stats"
    }

    output {
        File mergedStats = "merged.stats"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        stats: {description: "Statistics files to be merged.", category: "required"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        mergedStats: {description: "Merged stats from scattered Mutect2 runs."}
    }
}

task ModelSegments {
    input {
        String outputDir = "."
        String outputPrefix
        File denoisedCopyRatios
        File allelicCounts
        Int minimumTotalAlleleCountCase = if defined(normalAllelicCounts) then 0 else 30
        Int maximumNumberOfSmoothingIterations = 10

        File? normalAllelicCounts

        String javaXmx = "10G"
        String memory = "11GiB"
        Int timeMinutes = 60
        String dockerImage = "quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0"
    }

    command {
        set -e
        mkdir -p ~{outputDir}
        gatk --java-options '-Xmx~{javaXmx} -XX:ParallelGCThreads=1' \
        ModelSegments \
        --denoised-copy-ratios ~{denoisedCopyRatios} \
        --allelic-counts ~{allelicCounts} \
        ~{"--normal-allelic-counts " + normalAllelicCounts} \
        --minimum-total-allele-count-case ~{minimumTotalAlleleCountCase} \
        --maximum-number-of-smoothing-iterations ~{maximumNumberOfSmoothingIterations} \
        --output ~{outputDir} \
        --output-prefix ~{outputPrefix}
    }

    output {
        File hetrozygousAllelicCounts = outputDir + "/" + outputPrefix + ".hets.tsv"
        File copyRatioSegments = outputDir + "/" + outputPrefix + ".cr.seg"
        File copyRatioCBS = outputDir + "/" + outputPrefix + ".cr.igv.seg"
        File alleleFractionCBS = outputDir + "/" + outputPrefix + ".af.igv.seg"
        File unsmoothedModeledSegments = outputDir + "/" + outputPrefix + ".modelBegin.seg"
        File unsmoothedCopyRatioParameters = outputDir + "/" + outputPrefix + ".modelBegin.cr.param"
        File unsmoothedAlleleFractionParameters = outputDir + "/" + outputPrefix + ".modelBegin.af.param"
        File modeledSegments = outputDir + "/" + outputPrefix + ".modelFinal.seg"
        File copyRatioParameters = outputDir + "/" + outputPrefix + ".modelFinal.cr.param"
        File alleleFractionParameters = outputDir + "/" + outputPrefix + ".modelFinal.af.param"
        File? normalHetrozygousAllelicCounts = outputDir + "/" + outputPrefix + ".hets.normal.tsv"
    }

    runtime {
        memory: memory
        time_minute: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        outputDir: {description: "The directory to write the ouput to.", category: "common"}
        outputPrefix: {description: "The prefix of the output files. Should not include directories.", category: "required"}
        denoisedCopyRatios: {description: "The denoised copy ratios as generated by DenoiseReadCounts.", category: "required"}
        allelicCounts: {description: "The allelicCounts as generate by CollectAllelicCounts.", category: "required" }
        minimumTotalAlleleCountCase: {description: "Equivalent to gatk ModelSeqments' `--minimum-total-allele-count-case` option.", category: "advanced"}
        maximumNumberOfSmoothingIterations: {description: "Equivalent to gatk ModelSeqments' `--maximum-number-of-smoothing-iterations` option.", category: "advanced"}
        normalAllelicCounts: {description: "The allelicCounts as generate by CollectAllelicCounts for a matched normal.", category: "common"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        hetrozygousAllelicCounts: {description: "Allelic-counts file containing the counts at sites genotyped as heterozygous in the case sample."}
        copyRatioSegments: {description: "It contains the segments from the .modelFinal.seg file converted to a format suitable for input to CallCopyRatioSegments."}
        copyRatioCBS: {description: "The posterior medians of the log2 copy ratio."}
        alleleFractionCBS: {description: "Minor-allele fraction."}
        unsmoothedModeledSegments: {description: "The initial modeled-segments result before segmentation smoothing."}
        unsmoothedCopyRatioParameters: {description: "The initial copy-ratio-model global-parameter result before segmentation smoothing."}
        unsmoothedAlleleFractionParameters: {description: "The initial allele-fraction-model global-parameter result before segmentation smoothing."}
        modeledSegments: {description: "The final modeled-segments result after segmentation smoothing."}
        copyRatioParameters: {description: "The final copy-ratio-model global-parameter result after segmentation smoothing."}
        alleleFractionParameters: {description: "The final allele-fraction-model global-parameter result after segmentation smoothing."}
        normalHetrozygousAllelicCounts: {description: "Allelic-counts file containing the counts at sites genotyped as heterozygous in the matched-normal sample."}
    }
}

task MuTect2 {
    input {
        Array[File]+ inputBams
        Array[File]+ inputBamsIndex
        File referenceFasta
        File referenceFastaDict
        File referenceFastaFai
        String outputVcf
        String tumorSample
        String f1r2TarGz = "f1r2.tar.gz"
        Array[File]+ intervals
        String outputStats = outputVcf + ".stats"

        String? normalSample
        File? germlineResource
        File? germlineResourceIndex
        File? panelOfNormals
        File? panelOfNormalsIndex

        String javaXmx = "4G"
        String memory = "5GiB"
        Int timeMinutes = 240
        String dockerImage = "quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputVcf})"
        gatk --java-options '-Xmx~{javaXmx} -XX:ParallelGCThreads=1' \
        Mutect2 \
        -R ~{referenceFasta} \
        -I ~{sep=" -I " inputBams} \
        -tumor ~{tumorSample} \
        ~{"-normal " + normalSample} \
        ~{"--germline-resource " + germlineResource} \
        ~{"--panel-of-normals " + panelOfNormals} \
        ~{"--f1r2-tar-gz " + f1r2TarGz} \
        -O ~{outputVcf} \
        -L ~{sep=" -L " intervals}
    }

    output {
        File vcfFile = outputVcf
        File vcfFileIndex = outputVcf + ".tbi"
        File f1r2File = f1r2TarGz
        File stats = outputStats
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        inputBams: {description: "The BAM files on which to perform variant calling.", category: "required"}
        inputBamsIndex: {description: "The indexes for the input BAM files.", category: "required"}
        referenceFasta: {description: "The reference fasta file which was also used for mapping.", category: "required"}
        referenceFastaDict: {description: "The sequence dictionary associated with the reference fasta file.", category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        outputVcf: {description: "The location to write the output VCF file to.", category: "required"}
        tumorSample: {description: "The name of the tumor/case sample.", category: "required"}
        f1r2TarGz: {description: "Equivalent to Mutect2's `--f1r2-tar-gz` option.", category: "advanced"}
        intervals: {description: "Bed files describing the regiosn to operate on.", category: "required"}
        outputStats: {description: "The location the output statistics should be written to.", category: "advanced"}
        normalSample: {description: "The name of the normal/control sample.", category: "common"}
        germlineResource: {description: "Equivalent to Mutect2's `--germline-resource` option.", category: "advanced"}
        germlineResourceIndex: {description: "The index for the germline resource.", category: "advanced"}
        panelOfNormals: {description: "Equivalent to Mutect2's `--panel-of-normals` option.", category: "advanced"}
        panelOfNormalsIndex: {description: "The index for the panel of normals.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        vcfFile: {description: "Somatic SNVs and indels called via local assembly of haplotypes."}
        vcfFileIndex: {description: "Index for Mutect2 VCF."}
        f1r2File: {description: "Contains information that can then be passed to LearnReadOrientationModel, which generate an artifact prior table for each tumor sample for FilterMutectCalls to use."}
        stats: {description: "Stats file."}
    }
}

task PlotDenoisedCopyRatios {
    input {
        File referenceFastaDict
        String outputDir = "."
        String outputPrefix
        File standardizedCopyRatios
        File denoisedCopyRatios

        Int? minimumContigLength

        String javaXmx = "3G"
        String memory = "4GiB"
        Int timeMinutes = 2
        String dockerImage = "broadinstitute/gatk:4.1.8.0"
    }

    command {
        set -e
        mkdir -p ~{outputDir}
        gatk --java-options '-Xmx~{javaXmx} -XX:ParallelGCThreads=1' \
        PlotDenoisedCopyRatios \
        --standardized-copy-ratios ~{standardizedCopyRatios} \
        --denoised-copy-ratios ~{denoisedCopyRatios} \
        --sequence-dictionary ~{referenceFastaDict} \
        ~{"--minimum-contig-length " + minimumContigLength} \
        --output ~{outputDir} \
        --output-prefix ~{outputPrefix}
    }

    output {
        File denoisedCopyRatiosPlot = outputDir + "/" + outputPrefix + ".denoised.png"
        File standardizedMedianAbsoluteDeviation = outputDir + "/" + outputPrefix + ".standardizedMAD.txt"
        File denoisedMedianAbsoluteDeviation = outputDir + "/" + outputPrefix + ".denoisedMAD.txt"
        File deltaMedianAbsoluteDeviation = outputDir + "/" + outputPrefix + ".deltaMAD.txt"
        File deltaScaledMedianAbsoluteDeviation = outputDir + "/" + outputPrefix + ".scaledDeltaMAD.txt"
        File? denoisedCopyRatiosLimitedPlot = outputDir + "/" + outputPrefix + ".denoisedLimit4.png"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        referenceFastaDict: {description: "The sequence dictionary associated with the reference fasta file used for the analyses.", category: "required"}
        outputDir: {description: "The directory to write the ouput to.", category: "common"}
        outputPrefix: {description: "The prefix of the output files. Should not include directories.", category: "required"}
        standardizedCopyRatios: {description: "The standardized copy ratios as generated by DenoiseReadCounts.", category: "required"}
        denoisedCopyRatios: {description: "The denoised copy ratios as generated by DenoiseReadCounts.", category: "required"}
        minimumContigLength: {description: "The minimum length for a contig to be included in the plots.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        denoisedCopyRatiosPlot: {description: "Plot showing the entire range of standardized and denoised copy ratios."}
        standardizedMedianAbsoluteDeviation: {description: "Standardized median absolute deviation copy ratios."}
        denoisedMedianAbsoluteDeviation: {description: "Denoised median absolute deviation copy ratios."}
        deltaMedianAbsoluteDeviation: {description: "The change between `standardizedMedianAbsoluteDeviation` & `denoisedMedianAbsoluteDeviation`."}
        deltaScaledMedianAbsoluteDeviation: {description: "The change between `standardizedMedianAbsoluteDeviation` & `denoisedMedianAbsoluteDeviation` scaled by standardized MAD."}
        denoisedCopyRatiosLimitedPlot: {description: "Plot showing the standardized and denoised copy ratios limited to ratios within [0, 4]."}
    }
}

task PlotModeledSegments {
    input {
        File referenceFastaDict
        String outputDir = "."
        String outputPrefix
        File denoisedCopyRatios
        File segments
        File allelicCounts

        Int? minimumContigLength

        String javaXmx = "3G"
        String memory = "4GiB"
        Int timeMinutes = 2
        String dockerImage = "broadinstitute/gatk:4.1.8.0"
    }

    command {
        set -e
        mkdir -p ~{outputDir}
        gatk --java-options '-Xmx~{javaXmx} -XX:ParallelGCThreads=1' \
        PlotModeledSegments \
        --denoised-copy-ratios ~{denoisedCopyRatios} \
        --allelic-counts ~{allelicCounts} \
        --segments ~{segments} \
        --sequence-dictionary ~{referenceFastaDict} \
        ~{"--minimum-contig-length " + minimumContigLength} \
        --output ~{outputDir} \
        --output-prefix ~{outputPrefix}
    }

    output {
        File modeledSegmentsPlot = outputDir + "/" + outputPrefix + ".modeled.png"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        referenceFastaDict: {description: "The sequence dictionary associated with the reference fasta file used for the analyses.", category: "required"}
        outputDir: {description: "The directory to write the ouput to.", category: "common"}
        outputPrefix: {description: "The prefix of the output files. Should not include directories.", category: "required"}
        denoisedCopyRatios: {description: "The denoised copy ratios as generated by DenoiseReadCounts.", category: "required"}
        segments: {description: "The modeled segments as generated by ModelSegments.", category: "required"}
        allelicCounts: {description: "The hetrozygous allelic counts as generated by ModelSegments.", category: "required"}
        minimumContigLength: {description: "The minimum length for a contig to be included in the plots.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        modeledSegmentsPlot: {description: "This plot shows the input denoised copy ratios and/or alternate-allele fractions as points, as well as box plots for the available posteriors in each segment."}
    }
}

task PreprocessIntervals {
    input {
        File referenceFasta
        File referenceFastaDict
        File referenceFastaFai
        String outputIntervalList = "bins.interval_list"
        Int binLength = if defined(intervals) then 0 else 1000
        Int padding = if defined(intervals) then 250 else 0
        String intervalMergingRule = "OVERLAPPING_ONLY"

        File? intervals

        String javaXmx = "3G"
        String memory = "4GiB"
        Int timeMinutes = 1 + ceil(size(referenceFasta, "G") * 6)
        String dockerImage = "quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputIntervalList})"
        gatk --java-options '-Xmx~{javaXmx} -XX:ParallelGCThreads=1' \
        PreprocessIntervals \
        -R ~{referenceFasta} \
        --sequence-dictionary ~{referenceFastaDict} \
        --bin-length ~{binLength} \
        --padding ~{padding} \
        ~{"-L " + intervals} \
        --interval-merging-rule ~{intervalMergingRule} \
        -O ~{outputIntervalList}
    }

    output {
        File intervalList = outputIntervalList
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        referenceFasta: {description: "The reference fasta file.", category: "required"}
        referenceFastaDict: {description: "The sequence dictionary associated with the reference fasta file.", category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        outputIntervalList: {description: "The location the output should be written to.", category: "advanced"}
        binLength: {description: "The size of the bins to be created. Should be 0 for targeted/exome sequencing.", category: "advanced"}
        padding: {description: "The padding to be added to the bins. Should be 0 if contiguos binning is used, eg with WGS.", category: "advanced"}
        intervalMergingRule: {description: "Equivalent to gatk PreprocessIntervals' `--interval-merging-rule` option.", category: "advanced"}
        intervals: {description: "Bed files describing the regiosn to operate on.", category: "common"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        intervalList: {description: "Preprocessed Picard interval-list file."}
    }
}

task SelectVariants {
    input {
        File inputVcf
        File inputVcfIndex
        File referenceFasta
        File referenceFastaDict
        File referenceFastaFai
        String outputPath = "output.vcf.gz"
        Array[File] intervals = []

        String? selectTypeToInclude

        String javaXmx = "4G"
        String memory = "5GiB"
        Int timeMinutes = 60
        String dockerImage = "quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPath})"
        gatk --java-options '-Xmx~{javaXmx} -XX:ParallelGCThreads=1' \
        SelectVariants \
        -R ~{referenceFasta} \
        -V ~{inputVcf} \
        ~{"--select-type-to-include " + selectTypeToInclude} \
        ~{true="-L" false="" length(intervals) > 0} ~{sep=' -L ' intervals} \
        -O ~{outputPath}
    }

    output {
        File outputVcf = outputPath
        File outputVcfIndex = outputPath + ".tbi"
    }

    runtime {
        memory: memory
        time_minute: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        inputVcf: {description: "The VCF input file.", category: "required"}
        inputVcfIndex: {description: "The input VCF file's index.", category: "required"}
        referenceFasta: {description: "The reference fasta file which was also used for mapping.", category: "required"}
        referenceFastaDict: {description: "The sequence dictionary associated with the reference fasta file.", category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        outputPath: {description: "The location the output VCF file should be written.", category: "advanced"}
        intervals: {description: "Bed files or interval lists describing the regions to operate on.", category: "common"}
        selectTypeToInclude: {description: "Select only a certain type of variants from the input file.", category: "common"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        outputVcf: {description: "A new VCF file containing the selected subset of variants."}
        outputVcfIndex: {description: "Index of the new output VCF file."}
    }
}

task SplitNCigarReads {
    input {
        File inputBam
        File inputBamIndex
        File referenceFasta
        File referenceFastaDict
        File referenceFastaFai
        String outputBam
        Array[File] intervals = []

        String javaXmx = "4G"
        String memory = "5GiB"
        Int timeMinutes = 120 # This will likely be used with intervals, as such size based estimation can't be used.
        String dockerImage = "quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputBam})"
        gatk --java-options '-Xmx~{javaXmx} -XX:ParallelGCThreads=1' \
        SplitNCigarReads \
        -I ~{inputBam} \
        -R ~{referenceFasta} \
        -O ~{outputBam} \
        ~{true="-L" false="" length(intervals) > 0} ~{sep=' -L ' intervals}
    }

    output {
        File bam = outputBam
        File bamIndex = sub(outputBam, "\.bam$", ".bai")
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        inputBam: {description: "The BAM file for which spliced reads should be split.", category: "required"}
        inputBamIndex: {description: "The input BAM file's index.", category: "required"}
        referenceFasta: {description: "The reference fasta file which was also used for mapping.", category: "required"}
        referenceFastaDict: {description: "The sequence dictionary associated with the reference fasta file.", category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        outputBam: {description: "The location the output BAM file should be written.", category: "required"}
        intervals: {description: "Bed files or interval lists describing the regions to operate on.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        bam: {description: "BAM file with reads split at N CIGAR elements and CIGAR strings updated."}
        bamIndex: {description: "Index of output BAM file."}
    }
}

task VariantEval {
    input {
        Array[File] evalVcfs
        Array[File] evalVcfsIndex
        Array[File] comparisonVcfs = []
        Array[File] comparisonVcfsIndex = []
        Array[File] intervals = []
        String outputPath = "eval.table"
        Boolean doNotUseAllStandardModules = false 
        Boolean doNotUseAllStandardStratifications = false 
        Array[String] evalModules = []
        Array[String] stratificationModules = []
        Array[String] samples = []
        Boolean mergeEvals = false

        File? referenceFasta
        File? referenceFastaDict
        File? referenceFastaFai
        File? dbsnpVCF
        File? dbsnpVCFIndex

        String javaXmx = "4G"
        String memory = "5GiB"
        # TODO: Refine estimate. For now 4 minutes per GB of input.
        Int timeMinutes = ceil(size(flatten([evalVcfs, comparisonVcfs, select_all([referenceFasta, dbsnpVCF])]), "G") * 20)
        String dockerImage = "quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPath})"
        gatk --java-options '-Xmx~{javaXmx} -XX:ParallelGCThreads=1' \
        VariantEval \
        --output ~{outputPath} \
        ~{true="--eval" false="" length(evalVcfs) > 0} ~{sep=" --eval " evalVcfs} \
        ~{true="--comparison" false="" length(comparisonVcfs) > 0} ~{sep=" --comparison " comparisonVcfs} \
        ~{"-R " + referenceFasta} \
        ~{"--dbsnp " + dbsnpVCF } \
        ~{true="-L" false="" length(intervals) > 0} ~{sep=' -L ' intervals} \
        ~{true="--sample" false="" length(samples) > 0} ~{sep=' --sample ' samples} \
        ~{true="--do-not-use-all-standard-modules" false="" doNotUseAllStandardModules} \
        ~{true="--do-not-use-all-standard-stratifications" false="" doNotUseAllStandardStratifications} \
        ~{true="-EV" false="" length(evalModules) > 0} ~{sep=" -EV " evalModules} \
        ~{true="-ST" false="" length(stratificationModules) > 0} ~{sep=" -ST " stratificationModules} \
        ~{true="--merge-evals" false="" mergeEvals}
    }

    output {
        File table = outputPath
    }

    runtime {
        cpu: 1
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        evalVcfs: {description: "Variant sets to evaluate.", category: "required"}
        evalVcfsIndex: {description: "Indexes for the variant sets.", category: "required"}
        comparisonVcfs: {description: "Compare set vcfs.", category: "advanced"}
        comparisonVcfsIndex: {description: "Indexes for the compare sets.", category: "advanced"}
        intervals: {description: "Bed files or interval lists describing the regions to operate on.", category: "advanced"}
        outputPath: {description: "The location the output table should be written.", category: "advanced"}
        doNotUseAllStandardModules: {description: "Do not use the standard modules by default (instead, only those that are specified with the evalModules option).", category: "common"}
        doNotUseAllStandardStratifications: {description: "Do not use the standard stratification modules by default (instead, only those that are specified with the stratificationModules option).", category: "common"}
        evalModules: {description: "One or more specific eval modules to apply to the eval track(s) (in addition to the standard modules, unless doNotUseAllStandardModules=true).", category: "common"}
        stratificationModules: {description: "One or more specific stratification modules to apply to the eval track(s) (in addition to the standard stratifications, unless doNotUseAllStandardStratifications=true).", category: "common"}
        samples: {description: "Derive eval and comp contexts using only these sample genotypes, when genotypes are available in the original context." , category: "advanced"}
        mergeEvals: {description: "If provided, all evalVcf tracks will be merged into a single eval track.", category: "common"}
        referenceFasta: {description: "The reference fasta file which was also used for mapping.", category: "common"}
        referenceFastaDict: {description: "The sequence dictionary associated with the reference fasta file.", category: "common"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "common"}
        dbsnpVCF: {description: "A dbSNP VCF.", category: "common"}
        dbsnpVCFIndex: {description: "The index for the dbSNP VCF.", category: "common"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        table: {description: "Evaluation tables detailing the results of the eval modules which were applied."}
    }
}

task VariantFiltration {
    input {
        File inputVcf
        File inputVcfIndex
        File referenceFasta
        File referenceFastaDict
        File referenceFastaFai
        String outputPath = "filtered.vcf.gz"
        Array[String]+ filterArguments
        Array[File] intervals = []

        String javaXmx = "4G"
        String memory = "5GiB"
        Int timeMinutes = 120
        String dockerImage = "quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPath})"
        gatk --java-options '-Xmx~{javaXmx} -XX:ParallelGCThreads=1' \
        VariantFiltration \
        -I ~{inputVcf} \
        -R ~{referenceFasta} \
        -O ~{outputPath} \
        ~{sep=" " filterArguments} \
        ~{true="-L" false="" length(intervals) > 0} ~{sep=' -L ' intervals}
    }

    output {
        File filteredVcf = outputPath
        File filteredVcfIndex = outputPath + ".tbi"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        inputVcf: {description: "The VCF to be filtered.", category: "required"}
        inputVcfIndex: {description: "The input VCF file's index.", category: "required"}
        referenceFasta: {description: "The reference fasta file which was also used for mapping.", category: "required"}
        referenceFastaDict: {description: "The sequence dictionary associated with the reference fasta file.", category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        outputPath: {description: "The location the output VCF file should be written.", category: "common"}
        filterArguments: {description: "Arguments that should be used for the filter. For example: ['--filter-name', 'my_filter', '--filter-expression', 'AB<0.2'].", category: "required"}
        intervals: {description: "Bed files or interval lists describing the regions to operate on.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        filteredVcf: {description: "A filtered VCF in which passing variants are annotated as PASS and failing variants are annotated with the name(s) of the filter(s) they failed."}
        filteredVcfIndex: {description: "Index of filtered VCF."}
    }
}
