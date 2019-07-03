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
        File? panelOfNormals
        Array[File]+ intervals

        Int memory = 4
        Float memoryMultiplier = 3
        String dockerImage = "quay.io/biocontainers/gatk4:4.1.0.0--0"
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
        ~{"--panel-of-normals " + panelOfNormals} \
        -O ~{outputVcf} \
        -L ~{sep=" -L " intervals}
    }

    output {
        File vcfFile = outputVcf
        File vcfFileIndex = outputVcf + ".tbi"
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
