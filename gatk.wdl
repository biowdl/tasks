version 1.0

# Apply Base Quality Score Recalibration (BQSR) model
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

        String memory = "12G"
        String javaXmx = "4G"
        String dockerImage = "quay.io/biocontainers/gatk4:4.1.0.0--0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputBamPath})"
        gatk --java-options -Xmx~{javaXmx} \
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
        docker: dockerImage
        memory: memory
    }

    parameter_meta {
        inputBam: {description: "The BAM file which should be recalibrated.", category: "required"}
        inputBamIndex: {description: "The input BAM file's index.", category: "required"}
        outputBamPath: {description: "The location the resulting BAM file should be written.", category: "required"}
        recalibrationReport: {description: "The BQSR report the be used for recalibration.", category: "required"}
        sequenceGroupInterval: {description: "Bed files describing the regions to operate on.", category: "advanced"}
        referenceFasta: {description: "The reference fasta file which was also used for mapping.",
                         category: "required"}
        referenceFastaDict: {description: "The sequence dictionary associated with the reference fasta file.",
                             category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}

        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.",
                  category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}

# Generate Base Quality Score Recalibration (BQSR) model
task BaseRecalibrator {
    input {
        File inputBam
        File inputBamIndex
        String recalibrationReportPath
        Array[File] sequenceGroupInterval = []
        Array[File] knownIndelsSitesVCFs = []
        Array[File] knownIndelsSitesVCFIndexes = []
        File? dbsnpVCF
        File? dbsnpVCFIndex
        File referenceFasta
        File referenceFastaDict
        File referenceFastaFai

        String memory = "12G"
        String javaXmx = "4G"
        String dockerImage = "quay.io/biocontainers/gatk4:4.1.0.0--0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{recalibrationReportPath})"
        gatk --java-options -Xmx~{javaXmx} \
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
        docker: dockerImage
        memory: memory
    }

    parameter_meta {
        inputBam: {description: "The BAM file to generate a BQSR report for.", category: "required"}
        inputBamIndex: {description: "The index of the input BAM file.", category: "required"}
        recalibrationReportPath: {description: "The location to write the BQSR report to.", category: "required"}
        sequenceGroupInterval: {description: "Bed files describing the regions to operate on.", category: "advanced"}
        knownIndelsSitesVCFs: {description: "VCF files with known indels.", category: "advanced"}
        knownIndelsSitesVCFIndexes: {description: "The indexed for the known variant VCFs.", category: "advanced"}
        dbsnpVCF: {description: "A dbSNP VCF.", category: "common"}
        dbsnpVCFIndex: {description: "The index for the dbSNP VCF.", category: "common"}
        referenceFasta: {description: "The reference fasta file which was also used for mapping.",
                         category: "required"}
        referenceFastaDict: {description: "The sequence dictionary associated with the reference fasta file.",
                             category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}

        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.",
                  category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
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

        String memory = "24G"
        String javaXmx = "12G"
        String dockerImage = "quay.io/biocontainers/gatk4:4.1.0.0--0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPath})"
        gatk --java-options -Xmx~{javaXmx} \
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
        docker: dockerImage
        memory: memory
    }

    parameter_meta {
        gvcfFiles: {description: "The GVCF files to be combined.", category: "required"}
        gvcfFilesIndex: {description: "The indexes for the GVCF files.", caregory: "required"}
        intervals: {description: "Bed files or interval lists describing the regions to operate on.", category: "advanced"}
        outputPath: {description: "The location the combined GVCF should be written to.", category: "required"}
        referenceFasta: {description: "The reference fasta file which was also used for mapping.",
                         category: "required"}
        referenceFastaDict: {description: "The sequence dictionary associated with the reference fasta file.",
                             category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}

        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.",
                  category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}

# Combine multiple recalibration tables from scattered BaseRecalibrator runs
task GatherBqsrReports {
    input {
        Array[File] inputBQSRreports
        String outputReportPath

        String memory = "12G"
        String javaXmx = "4G"
        String dockerImage = "quay.io/biocontainers/gatk4:4.1.0.0--0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputReportPath})"
        gatk --java-options -Xmx~{javaXmx} \
        GatherBQSRReports \
        -I ~{sep=' -I ' inputBQSRreports} \
        -O ~{outputReportPath}
    }

    output {
        File outputBQSRreport = outputReportPath
    }

    runtime {
        docker: dockerImage
        memory: memory
    }

    parameter_meta {
        inputBQSRreports: {description: "The BQSR reports to be merged.", category: "required"}
        outputReportPath: {description: "The location of the combined BQSR report.", category: "required"}

        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.",
                  category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}

task GenotypeGVCFs {
    input {
        Array[File]+ gvcfFiles
        Array[File]+ gvcfFilesIndex
        Array[File]+ intervals
        String outputPath
        File referenceFasta
        File referenceFastaDict
        File referenceFastaFai
        File? dbsnpVCF
        File? dbsnpVCFIndex

        String memory = "18G"
        String javaXmx = "6G"
        String dockerImage = "quay.io/biocontainers/gatk4:4.1.0.0--0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPath})"
        gatk --java-options -Xmx~{javaXmx} \
        GenotypeGVCFs \
        -R ~{referenceFasta} \
        -O ~{outputPath} \
        ~{true="-D" false="" defined(dbsnpVCF)} ~{dbsnpVCF} \
        -G StandardAnnotation \
        --only-output-calls-starting-in-intervals \
        -new-qual \
        -V ~{sep=' -V ' gvcfFiles} \
        -L ~{sep=' -L ' intervals}
    }

    output {
        File outputVCF = outputPath
        File outputVCFIndex = outputPath + ".tbi"

    }

    runtime {
        docker: dockerImage
        memory: memory
    }

    parameter_meta {
        gvcfFiles: {description: "The GVCF files to be genotypes.", category: "required"}
        gvcfFilesIndex: {description: "The index of the input GVCF files.", category: "required"}
        intervals: {description: "Bed files or interval lists describing the regions to operate on.", category: "required"}
        outputPath: {description: "The location to write the output VCF file to.", category: "required"}
        referenceFasta: {description: "The reference fasta file which was also used for mapping.",
                         category: "required"}
        referenceFastaDict: {description: "The sequence dictionary associated with the reference fasta file.",
                             category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        dbsnpVCF: {description: "A dbSNP VCF.", category: "common"}
        dbsnpVCFIndex: {description: "The index for the dbSNP VCF.", category: "common"}

        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.",
                  category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}

# Call variants on a single sample with HaplotypeCaller to produce a GVCF
task HaplotypeCallerGvcf {
    input {
        Array[File]+ inputBams
        Array[File]+ inputBamsIndex
        Array[File]+? intervalList
        Array[File]+? excludeIntervalList
        String gvcfPath
        File referenceFasta
        File referenceFastaIndex
        File referenceFastaDict
        Float contamination = 0.0
        File? dbsnpVCF
        File? dbsnpVCFIndex
        Int? ploidy

        String memory = "12G"
        String javaXmx = "4G"
        String dockerImage = "quay.io/biocontainers/gatk4:4.1.0.0--0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{gvcfPath})"
        gatk --java-options -Xmx~{javaXmx} \
        HaplotypeCaller \
        -R ~{referenceFasta} \
        -O ~{gvcfPath} \
        -I ~{sep=" -I " inputBams} \
        ~{"--sample-ploidy " + ploidy} \
        ~{true="-L" false="" defined(intervalList)} ~{sep=' -L ' intervalList} \
        ~{true="-XL" false="" defined(excludeIntervalList)} ~{sep=' -XL ' excludeIntervalList} \
        ~{true="-D" false="" defined(dbsnpVCF)} ~{dbsnpVCF} \
        -contamination ~{contamination} \
        -ERC GVCF
    }

    output {
        File outputGVCF = gvcfPath
        File outputGVCFIndex = gvcfPath + ".tbi"
    }

    runtime {
        docker: dockerImage
        memory: memory
    }

    parameter_meta {
        inputBams: {description: "The BAM files on which to perform variant calling.", category: "required"}
        inputBamsIndex: {description: "The indexes for the input BAM files.", category: "required"}
        intervalList: {description: "Bed files or interval lists describing the regions to operate on.", category: "required"}
        gvcfPath: {description: "The location to write the output GVCF to.", category: "required"}
        referenceFasta: {description: "The reference fasta file which was also used for mapping.",
                         category: "required"}
        referenceFastaDict: {description: "The sequence dictionary associated with the reference fasta file.",
                             category: "required"}
        referenceFastaIndex: {description: "The index for the reference fasta file.", category: "required"}
        contamination: {description: "Equivalent to HaplotypeCaller's `-contamination` option.", category: "advanced"}
        dbsnpVCF: {description: "A dbSNP VCF.", category: "common"}
        dbsnpVCFIndex: {description: "The index for the dbSNP VCF.", category: "common"}

        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.",
                  category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
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
        String? normalSample
        File? germlineResource
        File? germlineResourceIndex
        File? panelOfNormals
        File? panelOfNormalsIndex
        String f1r2TarGz = "f1r2.tar.gz"
        Array[File]+ intervals
        String outputStats = outputVcf + ".stats"

        String memory = "16G"
        String javaXmx = "4G"
        String dockerImage = "quay.io/biocontainers/gatk4:4.1.2.0--1"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputVcf})"
        gatk --java-options -Xmx~{javaXmx} \
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
        docker: dockerImage
        memory: memory
    }

    parameter_meta {
        inputBams: {description: "The BAM files on which to perform variant calling.", category: "required"}
        inputBamsIndex: {description: "The indexes for the input BAM files.", category: "required"}
        referenceFasta: {description: "The reference fasta file which was also used for mapping.", category: "required"}
        referenceFastaDict: {description: "The sequence dictionary associated with the reference fasta file.", category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        outputVcf: {description: "The location to write the output VCF file to.", category: "required"}
        tumorSample: {description: "The name of the tumor/case sample.", category: "required"}
        normalSample: {description: "The name of the normal/control sample.", category: "common"}
        germlineResource: {description: "Equivalent to Mutect2's `--germline-resource` option.", category: "advanced"}
        germlineResourceIndex: {description: "The index for the germline resource.", category: "advanced"}
        panelOfNormals: {description: "Equivalent to Mutect2's `--panel-of-normals` option.", category: "advanced"}
        panelOfNormalsIndex: {description: "The index for the panel of normals.", category: "advanced"}
        f1r2TarGz: {description: "Equivalent to Mutect2's `--f1r2-tar-gz` option.", category: "advanced"}
        intervals: {description: "Bed files describing the regiosn to operate on.", category: "required"}
        outputStats: {description: "The location the output statistics should be written to.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.",
                  category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}

task LearnReadOrientationModel {
    input {
        Array[File]+ f1r2TarGz

        String memory = "24G"
        String javaXmx = "12G"
        String dockerImage = "quay.io/biocontainers/gatk4:4.1.2.0--1"
    }

    command {
        set -e
        gatk --java-options -Xmx~{javaXmx} \
        LearnReadOrientationModel \
        -I ~{sep=" -I " f1r2TarGz} \
        -O "artifact-priors.tar.gz"
    }

    output {
        File artifactPriorsTable = "artifact-priors.tar.gz"
    }

    runtime {
        docker: dockerImage
        memory: memory
    }

    parameter_meta {
        f1r2TarGz: {description: "A f1r2TarGz file outputed by mutect2.", category: "required"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.",
                  category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}

task MergeStats {
    input {
        Array[File]+ stats

        String memory = "28G"
        String javaXmx = "14G"
        String dockerImage = "quay.io/biocontainers/gatk4:4.1.2.0--1"
    }

    command {
        set -e
        gatk --java-options -Xmx~{javaXmx} \
        MergeMutectStats \
        -stats ~{sep=" -stats " stats} \
        -O "merged.stats"
    }

    output {
        File mergedStats = "merged.stats"
    }

    runtime {
        docker: dockerImage
        memory: memory
    }

    parameter_meta {
        stats: {description: "Statistics files to be merged.", category: "required"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.",
                  category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
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

        String memory = "24G"
        String javaXmx = "12G"
        String dockerImage = "quay.io/biocontainers/gatk4:4.1.2.0--1"
    }

    command {
        set -e
        gatk --java-options -Xmx~{javaXmx} \
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
        docker: dockerImage
        memory: memory
    }

    parameter_meta {
        sampleBam: {description: "A BAM file for which a pileup should be created.", category: "required"}
        sampleBamIndex: {description: "The index of the input BAM file.", category: "required"}
        variantsForContamination: {description: "A VCF file with common variants.", category: "required"}
        variantsForContaminationIndex: {description: "The index for the common variants VCF file.", category: "required"}
        sitesForContamination: {description: "A bed file describing regions to operate on.", category: "required"}
        sitesForContaminationIndex: {description: "The index for the bed file.", category: "required"}
        outputPrefix: {description: "The prefix for the ouput.", category: "required"}

        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.",
                  category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}

task CalculateContamination {
    input {
        File tumorPileups
        File? normalPileups

        String memory = "24G"
        String javaXmx = "12G"
        String dockerImage = "quay.io/biocontainers/gatk4:4.1.2.0--1"
    }

    command {
        set -e
        gatk --java-options -Xmx~{javaXmx} \
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
        docker: dockerImage
        memory: memory
    }

    parameter_meta {
        tumorPileups: {description: "The pileup summary of a tumor/case sample.", category: "required"}
        normalPileups: {description: "The pileup summary of the normal/control sample.", category: "common"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.",
                  category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
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
        File? contaminationTable
        File? mafTumorSegments
        File? artifactPriors
        Int uniqueAltReadCount = 4
        File mutect2Stats

        String memory = "24G"
        String javaXmx = "12G"
        String dockerImage = "quay.io/biocontainers/gatk4:4.1.2.0--1"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputVcf})"
        gatk --java-options -Xmx~{javaXmx} \
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
        docker: dockerImage
        memory: memory
    }

    parameter_meta {
        referenceFasta: {description: "The reference fasta file which was also used for mapping.", category: "required"}
        referenceFastaDict: {description: "The sequence dictionary associated with the reference fasta file.", category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        unfilteredVcf: {description: "An unfiltered VCF file as produced by Mutect2.", category: "required"}
        unfilteredVcfIndex: {description: "The index of the unfiltered VCF file.", category: "required"}
        outputVcf: {description: "The location the filtered VCF file should be written.", category: "required"}
        contaminationTable: {description: "Equivalent to FilterMutectCalls' `--contamination-table` option.", category: "advanced"}
        mafTumorSegments: {description: "Equivalent to FilterMutectCalls' `--tumor-segmentation` option.", category: "advanced"}
        artifactPriors: {description: "Equivalent to FilterMutectCalls' `--ob-priors` option.", category: "advanced"}
        uniqueAltReadCount: {description: "Equivalent to FilterMutectCalls' `--unique-alt-read-count` option.", category: "advanced"}
        mutect2Stats: {description: "Equivalent to FilterMutectCalls' `-stats` option.", category: "advanced"}

        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.",
                  category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
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

        String memory = "16G"
        String javaXmx = "4G"
        String dockerImage = "quay.io/biocontainers/gatk4:4.1.0.0--0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputBam})"
        gatk --java-options -Xmx~{javaXmx} \
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
        docker: dockerImage
        memory: memory
    }

    parameter_meta {
        inputBam: {description: "The BAM file for which spliced reads should be split.", category: "required"}
        inputBamIndex: {description: "The input BAM file's index.", category: "required"}
        referenceFasta: {description: "The reference fasta file which was also used for mapping.",
                         category: "required"}
        referenceFastaDict: {description: "The sequence dictionary associated with the reference fasta file.",
                             category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        outputBam: {description: "The location the output BAM file should be written.", category: "required"}
        intervals: {description: "Bed files or interval lists describing the regions to operate on.", category: "advanced"}

        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.",
                  category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
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
        Array[File]+ variantVcfs # follow "identifiers" array order
        Array[File]+ variantIndexes
        String outputPath

        String memory = "24G"
        String javaXmx = "12G"
        String dockerImage = "broadinstitute/gatk3:3.8-1"
    }

    command <<<
        set -e
        mkdir -p "$(dirname ~{outputPath})"

        # build "-V:<ID> <file.vcf>" arguments according to IDs and VCFs to merge
        # Make sure commands are run in bash
        V_args=$(bash -c '
        set -eu
        ids=(~{sep=" " identifiers})
        vars=(~{sep=" " variantVcfs})
        for (( i = 0; i < ${#ids[@]}; ++i ))
          do
            printf -- "-V:%s %s " "${ids[i]}" "${vars[i]}"
          done
        ')
        java -Xmx~{javaXmx} -jar /usr/GenomeAnalysisTK.jar \
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
        docker: dockerImage
        memory: memory
    }

    parameter_meta {
        referenceFasta: {description: "The reference fasta file which was also used for mapping.", category: "required"}
        referenceFastaDict: {description: "The sequence dictionary associated with the reference fasta file.", category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        genotypeMergeOption: {description: "Equivalent to CombineVariants' `--genotypemergeoption` option.", category: "advanced"}
        filteredRecordsMergeType: {description: "Equivalent to CombineVariants' `--filteredrecordsmergetype` option.", category: "advanced"}
        identifiers: {description: "The sample identifiers in the same order as variantVcfs.", category: "required"}
        variantVcfs: {description: "The input VCF files in the same order as identifiers.", category: "required"}
        variantIndexes: {description: "The indexes of the input VCF files.", category: "required"}
        outputPath: {description: "The location the output should be written to", category: "required"}

        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.",
                  category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}
