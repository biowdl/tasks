version 1.0

# Apply Base Quality Score Recalibration (BQSR) model
task ApplyBQSR {
    input {
        File inputBam
        File inputBamIndex
        String outputBamPath
        File recalibrationReport
        Array[File]+ sequenceGroupInterval
        File referenceFasta
        File referenceFastaDict
        File referenceFastaFai

        Int memory = 4
        Float memoryMultiplier = 3.0
        String dockerImage = "quay.io/biocontainers/gatk4:4.1.0.0--0"
    }

    command {
        set -e
        mkdir -p $(dirname ~{outputBamPath})
        gatk --java-options -Xmx~{memory}G \
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
        -L ~{sep=" -L " sequenceGroupInterval}
    }

    output {
        File recalibratedBam = outputBamPath
        File recalibratedBamIndex = sub(outputBamPath, "\.bam$", ".bai")
        File recalibratedBamMd5 = outputBamPath + ".md5"
    }

    runtime {
        docker: dockerImage
        memory: ceil(memory * memoryMultiplier)
    }
}

# Generate Base Quality Score Recalibration (BQSR) model
task BaseRecalibrator {
    input {
        File inputBam
        File inputBamIndex
        String recalibrationReportPath
        Array[File]+ sequenceGroupInterval
        Array[File]? knownIndelsSitesVCFs
        Array[File]? knownIndelsSitesVCFIndexes
        File? dbsnpVCF
        File? dbsnpVCFIndex
        File referenceFasta
        File referenceFastaDict
        File referenceFastaFai

        Int memory = 4
        Float memoryMultiplier = 3.0
        String dockerImage = "quay.io/biocontainers/gatk4:4.1.0.0--0"
    }

    Array[File]+ knownIndelsSitesVCFsArg = flatten([
        select_first([knownIndelsSitesVCFs, []]),
        [select_first([dbsnpVCF])]
    ])

    command {
        set -e
        mkdir -p $(dirname ~{recalibrationReportPath})
        gatk --java-options -Xmx~{memory}G \
        BaseRecalibrator \
        -R ~{referenceFasta} \
        -I ~{inputBam} \
        --use-original-qualities \
        -O ~{recalibrationReportPath} \
        --known-sites ~{sep=" --known-sites " knownIndelsSitesVCFsArg} \
        -L ~{sep=" -L " sequenceGroupInterval}
    }

    output {
        File recalibrationReport = recalibrationReportPath
    }

    runtime {
        docker: dockerImage
        memory: ceil(memory * memoryMultiplier)
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

        Int memory = 4
        Float memoryMultiplier = 3.0
        String dockerImage = "quay.io/biocontainers/gatk4:4.1.0.0--0"
    }

    command {
        set -e
        mkdir -p $(dirname ~{outputPath})
        gatk --java-options -Xmx~{memory}G \
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
        memory: ceil(memory * memoryMultiplier)
    }
}

# Combine multiple recalibration tables from scattered BaseRecalibrator runs
task GatherBqsrReports {
    input {
        Array[File] inputBQSRreports
        String outputReportPath

        Int memory = 4
        Float memoryMultiplier = 3.0
        String dockerImage = "quay.io/biocontainers/gatk4:4.1.0.0--0"
    }

    command {
        set -e
        mkdir -p $(dirname ~{outputReportPath})
        gatk --java-options -Xmx~{memory}G \
        GatherBQSRReports \
        -I ~{sep=' -I ' inputBQSRreports} \
        -O ~{outputReportPath}
    }

    output {
        File outputBQSRreport = outputReportPath
    }

    runtime {
        docker: dockerImage
        memory: ceil(memory * memoryMultiplier)
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
        Int memory = 6
        Float memoryMultiplier = 2.0
        String dockerImage = "quay.io/biocontainers/gatk4:4.1.0.0--0"
    }

    command {
        set -e
        mkdir -p $(dirname ~{outputPath})
        gatk --java-options -Xmx~{memory}G \
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
        memory: ceil(memory * memoryMultiplier)
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
        Int memory = 4
        Float memoryMultiplier = 3
        String dockerImage = "quay.io/biocontainers/gatk4:4.1.0.0--0"
    }

    command {
        set -e
        mkdir -p $(dirname ~{gvcfPath})
        gatk --java-options -Xmx~{memory}G \
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
        memory: ceil(memory * memoryMultiplier)
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

        Int memory = 4
        Float memoryMultiplier = 3
        String dockerImage = "quay.io/biocontainers/gatk4:4.1.2.0--1"
    }

    command {
        set -e
        mkdir -p $(dirname ~{outputVcf})
        gatk --java-options -Xmx~{memory}G \
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
        memory: ceil(memory * memoryMultiplier)
    }
}

task LearnReadOrientationModel {
    input {
        Array[File]+ f1r2TarGz

        Int memory = 8
        Float memoryMultiplier = 1.5
        String dockerImage = "quay.io/biocontainers/gatk4:4.1.2.0--1"
    }

    command {
        set -e
        gatk --java-options -Xmx~{memory}G \
        LearnReadOrientationModel \
        -I ~{sep=" -I " f1r2TarGz} \
        -O "artifact-priors.tar.gz"
    }

    output {
        File artifactPriorsTable = "artifact-priors.tar.gz"
    }

    runtime {
        docker: dockerImage
        memory: ceil(memory * memoryMultiplier)
    }
}

task MergeStats {
    input {
        Array[File]+ stats

        Int memory = 2
        Float memoryMultiplier = 1.5
        String dockerImage = "quay.io/biocontainers/gatk4:4.1.2.0--1"
    }

    command {
        set -e
        gatk --java-options -Xmx~{memory}G \
        MergeMutectStats \
        -stats ~{sep=" -stats " stats} \
        -O "merged.stats"
    }

    output {
        File mergedStats = "merged.stats"
    }

    runtime {
        docker: dockerImage
        memory: ceil(memory * memoryMultiplier)
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

        Int memory = 4
        Float memoryMultiplier = 1.5
        String dockerImage = "quay.io/biocontainers/gatk4:4.1.2.0--1"
    }

    command {
        set -e
        gatk --java-options -Xmx~{memory}G \
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
        memory: ceil(memory * memoryMultiplier)
    }
}

task CalculateContamination {
    input {
        File tumorPileups
        File? normalPileups

        Int memory = 4
        Float memoryMultiplier = 1.5
        String dockerImage = "quay.io/biocontainers/gatk4:4.1.2.0--1"
    }

    command {
        set -e
        gatk --java-options -Xmx~{memory}G \
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
        memory: ceil(memory * memoryMultiplier)
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

        Int memory = 4
        Float memoryMultiplier = 1.5
        String dockerImage = "quay.io/biocontainers/gatk4:4.1.2.0--1"
    }

    command {
        set -e
        mkdir -p $(dirname ~{outputVcf})
        gatk --java-options -Xmx~{memory}G \
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
        memory: ceil(memory * memoryMultiplier)
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
        Array[File]+ intervals

        Int memory = 4
        Float memoryMultiplier = 4
        String dockerImage = "quay.io/biocontainers/gatk4:4.1.0.0--0"
    }

    command {
        set -e
        mkdir -p $(dirname ~{outputBam})
        gatk --java-options -Xmx~{memory}G \
        SplitNCigarReads \
        -I ~{inputBam} \
        -R ~{referenceFasta} \
        -O ~{outputBam} \
        -L ~{sep=' -L ' intervals}
    }

    output {
        File bam = outputBam
        File bamIndex = sub(outputBam, "\.bam$", ".bai")
    }

    runtime {
        docker: dockerImage
        memory: ceil(memory * memoryMultiplier)
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

        Int memory = 4
        Float memoryMultiplier = 1.5
        String dockerImage = "broadinstitute/gatk3:3.8-1"
    }

    command <<<
        set -e -o pipefail
        mkdir -p $(dirname "~{outputPath}")

        # build "-V:<ID> <file.vcf>" arguments according to IDs and VCFs to merge
        ids=(~{sep=" " identifiers})
        vars=(~{sep=" " variantVcfs})
        V_args=$(
            for (( i = 0; i < ${#ids[@]}; ++i ))
              do
                printf -- "-V:%s %s " "${ids[i]}" "${vars[i]}"
              done
        )

        java -Xmx~{memory}G -jar ~{installDir}/GenomeAnalysisTK.jar \
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
        memory: ceil(memory * memoryMultiplier)
    }
}