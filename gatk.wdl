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
        mkdir -p $(dirname ~{outputBamPath})
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
        mkdir -p $(dirname ~{recalibrationReportPath})
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
}

task CombineGVCFs {
    input {
        Array[File]+ gvcfFiles
        Array[File]+ gvcfFilesIndex
        Array[File]+ intervals
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
        mkdir -p $(dirname ~{outputPath})
        gatk --java-options -Xmx~{javaXmx} \
        CombineGVCFs \
        -R ~{referenceFasta} \
        -O ~{outputPath} \
        -V ~{sep=' -V ' gvcfFiles} \
        -L ~{sep=' -L ' intervals}
    }

    output {
        File outputVcf = outputPath
        File outputVcfIndex = outputPath + ".tbi"
    }

    runtime {
        docker: dockerImage
        memory: memory
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
        mkdir -p $(dirname ~{outputReportPath})
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
        mkdir -p $(dirname ~{outputPath})
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
}

# Call variants on a single sample with HaplotypeCaller to produce a GVCF
task HaplotypeCallerGvcf {
    input {
        Array[File]+ inputBams
        Array[File]+ inputBamsIndex
        Array[File]+ intervalList
        String gvcfPath
        File referenceFasta
        File referenceFastaIndex
        File referenceFastaDict
        Float contamination = 0.0
        File? dbsnpVCF
        File? dbsnpVCFIndex

        String memory = "12G"
        String javaXmx = "4G"
        String dockerImage = "quay.io/biocontainers/gatk4:4.1.0.0--0"
    }

    command {
        set -e
        mkdir -p $(dirname ~{gvcfPath})
        gatk --java-options -Xmx~{javaXmx} \
        HaplotypeCaller \
        -R ~{referenceFasta} \
        -O ~{gvcfPath} \
        -I ~{sep=" -I " inputBams} \
        -L ~{sep=' -L ' intervalList} \
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
        mkdir -p $(dirname ~{outputVcf})
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
        String? extraArgs

        String memory = "24G"
        String javaXmx = "12G"
        String dockerImage = "quay.io/biocontainers/gatk4:4.1.2.0--1"
    }

    command {
        set -e
        mkdir -p $(dirname ~{outputVcf})
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
        --showHidden \
        ~{extraArgs}
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
        mkdir -p $(dirname ~{outputBam})
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
}

task CombineVariants {
    input {
        String installDir = "/usr"  # .jar location in the docker image

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
        mkdir -p $(dirname "~{outputPath}")

        # build "-V:<ID> <file.vcf>" arguments according to IDs and VCFs to merge
        # Make sure commands are run in bash
        bash -c '#!/usr/bin/env bash
        set -eux
        ids=(~{sep=" " identifiers})
        vars=(~{sep=" " variantVcfs})
        V_args=$(
            for (( i = 0; i < ${#ids[@]}; ++i ))
              do
                printf -- "-V:%s %s " "${ids[i]}" "${vars[i]}"
              done
        )
        java -Xmx~{javaXmx} -jar ~{installDir}/GenomeAnalysisTK.jar \
        -T CombineVariants \
        -R ~{referenceFasta} \
        --genotypemergeoption ~{genotypeMergeOption} \
        --filteredrecordsmergetype ~{filteredRecordsMergeType} \
        --out ~{outputPath} \
        $V_args
        '
    >>>

    output {
        File combinedVcf = outputPath
        File combinedVcfIndex = outputPath + ".tbi"
    }

    runtime {
        docker: dockerImage
        memory: memory
    }
}
