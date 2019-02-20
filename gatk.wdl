version 1.0

import "common.wdl"

# Apply Base Quality Score Recalibration (BQSR) model
task ApplyBQSR {
    input {
        IndexedBamFile inputBam
        String outputBamPath
        File recalibrationReport
        Array[File]+ sequenceGroupInterval
        Reference reference

        Int memory = 4
        Float memoryMultiplier = 3.0
        String dockerTag = "4.1.0.0--0"
    }

    command {
        set -e -o pipefail
        mkdir -p $(dirname ~{outputBamPath})
        gatk --java-options -Xmx~{memory}G \
        ApplyBQSR \
        --create-output-bam-md5 \
        --add-output-sam-program-record \
        -R ~{reference.fasta} \
        -I ~{inputBam.file} \
        --use-original-qualities \
        -O ~{outputBamPath} \
        -bqsr ~{recalibrationReport} \
        --static-quantized-quals 10 \
        --static-quantized-quals 20 \
        --static-quantized-quals 30 \
        -L ~{sep=" -L " sequenceGroupInterval}
    }

    output {
        IndexedBamFile recalibratedBam = {
            "file": outputBamPath,
            "index": sub(outputBamPath, "\.bam$", ".bai"),
            "md5": outputBamPath + ".md5"
        }
    }

    runtime {
        docker: "quay.io/biocontainers/gatk4:" + dockerTag
        memory: ceil(memory * memoryMultiplier)
    }
}

# Generate Base Quality Score Recalibration (BQSR) model
task BaseRecalibrator {
    input {
        IndexedBamFile inputBam
        String recalibrationReportPath
        Array[File]+ sequenceGroupInterval
        Array[File]? knownIndelsSitesVCFs
        Array[File]? knownIndelsSitesVCFIndexes
        IndexedVcfFile? dbsnpVCF
        Reference reference
        Int memory = 4
        Float memoryMultiplier = 3.0
        String dockerTag = "4.1.0.0--0"
    }

    Array[File]+ knownIndelsSitesVCFsArg = flatten([
        select_first([knownIndelsSitesVCFs, []]),
        [select_first([dbsnpVCF]).file]
    ])

    command {
        set -e -o pipefail
        mkdir -p $(dirname ~{recalibrationReportPath})
        gatk --java-options -Xmx~{memory}G \
        BaseRecalibrator \
        -R ~{reference.fasta} \
        -I ~{inputBam.file} \
        --use-original-qualities \
        -O ~{recalibrationReportPath} \
        --known-sites ~{sep=" --known-sites " knownIndelsSitesVCFsArg} \
        -L ~{sep=" -L " sequenceGroupInterval}
    }

    output {
        File recalibrationReport = recalibrationReportPath
    }

    runtime {
        docker: "quay.io/biocontainers/gatk4:" + dockerTag
        memory: ceil(memory * memoryMultiplier)
    }
}

task CombineGVCFs {
    input {
        Array[File]+ gvcfFiles
        Array[File]+ gvcfFilesIndex
        Array[File]+ intervals

        String outputPath


        Reference reference

        Int memory = 4
        Float memoryMultiplier = 3.0
        String dockerTag = "4.1.0.0--0"
    }

    command {
        set -e -o pipefail
        mkdir -p $(dirname ~{outputPath})
        gatk --java-options -Xmx~{memory}G \
        CombineGVCFs \
        -R ~{reference.fasta} \
        -O ~{outputPath} \
        -V ~{sep=' -V ' gvcfFiles} \
        -L ~{sep=' -L ' intervals}
    }

    output {
        IndexedVcfFile outputVCF = {
            "file": outputPath,
            "index": outputPath + ".tbi"
        }
    }

    runtime {
        docker: "quay.io/biocontainers/gatk4:" + dockerTag
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
        String dockerTag = "4.1.0.0--0"
    }

    command {
        set -e -o pipefail
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
        docker: "quay.io/biocontainers/gatk4:" + dockerTag
        memory: ceil(memory * memoryMultiplier)
    }
}

task GenotypeGVCFs {
    input {
        Array[File]+ gvcfFiles
        Array[File]+ gvcfFilesIndex
        Array[File]+ intervals

        String outputPath

        Reference reference

        IndexedVcfFile? dbsnpVCF

        Int memory = 6
        Float memoryMultiplier = 2.0
        String dockerTag = "4.1.0.0--0"
    }

    File dbsnpFile = if (defined(dbsnpVCF)) then select_first([dbsnpVCF]).file else ""

    command {
        set -e -o pipefail
        mkdir -p $(dirname ~{outputPath})
        gatk --java-options -Xmx~{memory}G \
        GenotypeGVCFs \
        -R ~{reference.fasta} \
        -O ~{outputPath} \
        ~{true="-D" false="" defined(dbsnpVCF)} ~{dbsnpFile} \
        -G StandardAnnotation \
        --only-output-calls-starting-in-intervals \
        -new-qual \
        -V ~{sep=' -V ' gvcfFiles} \
        -L ~{sep=' -L ' intervals}
    }

    output {
        IndexedVcfFile outputVCF = {
            "file": outputPath,
            "index": outputPath + ".tbi"
        }
    }

    runtime {
        docker: "quay.io/biocontainers/gatk4:" + dockerTag
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
        Reference reference
        Float contamination = 0.0

        IndexedVcfFile? dbsnpVCF

        Int memory = 4
        Float memoryMultiplier = 3
        String dockerTag = "4.1.0.0--0"
    }

    File dbsnpFile = if (defined(dbsnpVCF)) then select_first([dbsnpVCF]).file else ""

    command {
        set -e -o pipefail
        mkdir -p $(dirname ~{gvcfPath})
        gatk --java-options -Xmx~{memory}G \
        HaplotypeCaller \
        -R ~{reference.fasta} \
        -O ~{gvcfPath} \
        -I ~{sep=" -I " inputBams} \
        -L ~{sep=' -L ' intervalList} \
        ~{true="-D" false="" defined(dbsnpVCF)} ~{dbsnpFile} \
        -contamination ~{contamination} \
        -ERC GVCF
    }

    output {
        IndexedVcfFile outputGVCF = {
            "file": gvcfPath,
            "index": gvcfPath + ".tbi"
        }
    }

    runtime {
        docker: "quay.io/biocontainers/gatk4:" + dockerTag
        memory: ceil(memory * memoryMultiplier)
    }
}

task MuTect2 {
    input {
        Array[File]+ inputBams
        Array[File]+ inputBamsIndex
        Reference reference
        String outputVcf
        String tumorSample
        String? normalSample
        Array[File]+ intervals

        Int memory = 4
        Float memoryMultiplier = 3
        String dockerTag = "4.1.0.0--0"
    }

    command {
        set -e -o pipefail
        mkdir -p $(dirname ~{outputVcf})
        gatk --java-options -Xmx~{memory}G \
        Mutect2 \
        -R ~{reference.fasta} \
        -I ~{sep=" -I " inputBams} \
        -tumor ~{tumorSample} \
        ~{"-normal " + normalSample} \
        -O ~{outputVcf} \
        -L ~{sep=" -L " intervals}
    }

    output {
        IndexedVcfFile vcfFile = {
            "file": outputVcf,
            "index": outputVcf + ".tbi"
        }
    }

    runtime {
        docker: "quay.io/biocontainers/gatk4:" + dockerTag
        memory: ceil(memory * memoryMultiplier)
    }
}

task SplitNCigarReads {
    input {
        IndexedBamFile inputBam
        Reference reference
        String outputBam
        Array[File]+ intervals

        Int memory = 4
        Float memoryMultiplier = 4
        String dockerTag = "4.1.0.0--0"
    }

    command {
        set -e -o pipefail
        mkdir -p $(dirname ~{outputBam})
        gatk --java-options -Xmx~{memory}G \
        SplitNCigarReads \
        -I ~{inputBam.file} \
        -R ~{reference.fasta} \
        -O ~{outputBam} \
        -L ~{sep=' -L ' intervals}
    }

    output {
        IndexedBamFile bam = {
            "file": outputBam,
            "index": sub(outputBam, "\.bam$", ".bai")
        }
    }

    runtime {
        docker: "quay.io/biocontainers/gatk4:" + dockerTag
        memory: ceil(memory * memoryMultiplier)
    }
}
