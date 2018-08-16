version 1.0

# Apply Base Quality Score Recalibration (BQSR) model
task ApplyBQSR {
    input {
        String? preCommand
        File? gatkJar
        File inputBam
        File inputBamIndex
        String outputBamPath
        File recalibrationReport
        Array[File]+ sequenceGroupInterval
        File refDict
        File refFasta
        File refFastaIndex
        Int? compressionLevel

        Int memory = 4
        Float memoryMultiplier = 3.0
    }

    String toolCommand = if defined(gatkJar)
        then "java -Xmx" + memory + "G -jar " + gatkJar
        else "gatk --java-options -Xmx" + memory + "G"

    command {
        set -e -o pipefail
        ~{preCommand}
        ~{toolCommand} \
         ApplyBQSR \
         --create-output-bam-md5 \
         --add-output-sam-program-record \
         -R ~{refFasta} \
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
        File recalibrated_bam = outputBamPath
        File recalibrated_bam_checksum = outputBamPath + ".md5"
    }

    runtime {
        memory: ceil(memory * memoryMultiplier)
    }
}

# Generate Base Quality Score Recalibration (BQSR) model
task BaseRecalibrator {
    input {
        String? preCommand
        File? gatkJar
        File inputBam
        File inputBamIndex
        String recalibrationReportPath
        Array[File]+ sequenceGroupInterval
        Array[File]? knownIndelsSitesVCFs
        Array[File]? knownIndelsSitesIndices
        File? dbsnpVCF
        File? dbsnpVCFindex
        File refDict
        File refFasta
        File refFastaIndex
        Int memory = 4
        Float memoryMultiplier = 3.0
    }

    Array[File]+ knownIndelsSitesVCFsArg = flatten([
        select_first([knownIndelsSitesVCFs, []]),
        select_all([dbsnpVCF])
    ])
    Array[File]+ knownIndelsSitesIndicesArg = flatten([
        select_first([knownIndelsSitesIndices, []]),
        select_all([dbsnpVCFindex])
    ])

    String toolCommand = if defined(gatkJar)
        then "java -Xmx" + memory + "G -jar " + gatkJar
        else "gatk --java-options -Xmx" + memory + "G"

    command {
        set -e -o pipefail
        ~{preCommand}
        ~{toolCommand} \
        BaseRecalibrator \
        -R ~{refFasta} \
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
        memory: ceil(memory * memoryMultiplier)
    }
}

task CombineGVCFs {
    input {
        String? preCommand
        Array[File]+ gvcfFiles
        Array[File]+ gvcfFileIndexes
        Array[File]+ intervals

        String outputPath

        String? gatkJar

        File refFasta
        File refFastaIndex
        File refDict

        Int? compressionLevel #TODO This isn't being used?
        Int memory = 4
        Float memoryMultiplier = 3.0
    }

    String toolCommand = if defined(gatkJar)
        then "java -Xmx" + memory + "G -jar " + gatkJar
        else "gatk --java-options -Xmx" + memory + "G"

    command {
        set -e -o pipefail
        ~{preCommand}

        if [ ~{length(gvcfFiles)} -gt 1 ]; then
            ~{toolCommand} \
             CombineGVCFs \
             -R ~{refFasta} \
             -O ~{outputPath} \
             -V ~{sep=' -V ' gvcfFiles} \
             -L ~{sep=' -L ' intervals}
        else # TODO this should be handeled in wdl
            ln -sf ~{gvcfFiles[0]} ~{outputPath}
            ln -sf ~{gvcfFileIndexes[0]} ~{outputPath}.tbi
        fi
    }

    output {
        File outputGVCF = outputPath
        File outputGVCFindex = outputPath + ".tbi"
    }

    runtime {
        memory: ceil(memory * memoryMultiplier)
    }
}

# Combine multiple recalibration tables from scattered BaseRecalibrator runs
task GatherBqsrReports {
    input {
        String? preCommand
        String? gatkJar
        Array[File] inputBQSRreports
        String outputReportPath

        Int memory = 4
        Float memoryMultiplier = 3.0
    }

    String toolCommand = if defined(gatkJar)
        then "java -Xmx" + memory + "G -jar " + gatkJar
        else "gatk --java-options -Xmx" + memory + "G"

    command {
        set -e -o pipefail
        ~{preCommand}
        ~{toolCommand} \
        GatherBQSRReports \
        -I ~{sep=' -I ' inputBQSRreports} \
        -O ~{outputReportPath}
    }

    output {
        File outputBQSRreport = outputReportPath
    }

    runtime {
        memory: ceil(memory * memoryMultiplier)
    }
}

task GenotypeGVCFs {
    input {
        String? preCommand
        File gvcfFiles
        File gvcfFileIndexes
        Array[File]+ intervals

        String outputPath

        String? gatkJar

        File refFasta
        File refFastaIndex
        File refDict

        File? dbsnpVCF
        File? dbsnpVCFindex

        Int? compressionLevel
        Int memory = 4
        Float memoryMultiplier =3.0
    }

    String toolCommand = if defined(gatkJar)
        then "java -Xmx" + memory + "G -jar " + gatkJar
        else "gatk --java-options -Xmx" + memory + "G"

    command {
        set -e -o pipefail
        ~{preCommand}
        ~{toolCommand} \
        GenotypeGVCFs \
        -R ~{refFasta} \
        -O ~{outputPath} \
        ~{"-D " + dbsnpVCF} \
        -G StandardAnnotation \
        --only-output-calls-starting-in-intervals \
        -new-qual \
        -V ~{gvcfFiles} \
        -L ~{sep=' -L ' intervals}
    }

    output {
        File outputVCF = outputPath
        File outputVCFindex = outputPath + ".tbi"
    }

    runtime{
        memory: ceil(memory * memoryMultiplier)
    }
}

# Call variants on a single sample with HaplotypeCaller to produce a GVCF
task HaplotypeCallerGvcf {
    input {
        String? preCommand
        Array[File]+ inputBams
        Array[File]+ inputBamsIndex
        Array[File]+ intervalList
        String gvcfPath
        File refDict
        File refFasta
        File refFastaIndex
        Float contamination = 0.0
        Int? compressionLevel
        String? gatkJar

        File? dbsnpVCF
        File? dbsnpVCFindex

        Int memory = 4
        Float memoryMultiplier = 3
    }

    String toolCommand = if defined(gatkJar)
        then "java -Xmx" + memory + "G -jar " + gatkJar
        else "gatk --java-options -Xmx" + memory + "G"

    command {
        set -e -o pipefail
        ~{preCommand}
        ~{toolCommand} \
        HaplotypeCaller \
        -R ~{refFasta} \
        -O ~{gvcfPath} \
        -I ~{sep=" -I " inputBams} \
        -L ~{sep=' -L ' intervalList} \
        ~{"-D " + dbsnpVCF} \
        -contamination ~{contamination} \
        -ERC GVCF
    }

    output {
        File outputGVCF = gvcfPath
        File outputGVCFindex = gvcfPath + ".tbi"
    }

    runtime {
        memory: ceil(memory * memoryMultiplier)
    }
}

task MuTect2 {
    input {
        String? preCommand

        Array[File]+ inputBams
        Array[File]+ inputBamIndex
        File refFasta
        File refFastaIndex
        File refDict
        String outputVcf
        String tumorSample
        String? normalSample
        Array[File]+ intervals

        String? gatkJar
        Int memory = 4
        Float memoryMultiplier = 3
    }

    String toolCommand = if defined(gatkJar)
        then "java -Xmx" + memory + "G -jar " + gatkJar
        else "gatk --java-options -Xmx" + memory + "G"

    command {
        set -e -o pipefail
        ~{preCommand}
        ~{toolCommand} \
        Mutect2 \
        -R ~{refFasta} \
        -I ~{sep=" -I " inputBams} \
        -tumor ~{tumorSample} \
        ~{"-normal " + normalSample} \
        -O ~{outputVcf} \
        -L ~{sep=" -L " intervals}
    }

    output {
        File vcfFile = outputVcf
        File vcfIndex = outputVcf + ".tbi"
    }

    runtime {
        memory: ceil(memory * memoryMultiplier)
    }
}

task SplitNCigarReads {
    input {
        String? preCommand

        File inputBam
        File inputBamIndex
        File refFasta
        File refFastaIndex
        File refDict
        String outputBam
        String? gatkJar
        Array[File]+ intervals

        Int memory = 4
        Float memoryMultiplier = 3
    }

    String toolCommand = if defined(gatkJar)
        then "java -Xmx" + memory + "G -jar " + gatkJar
        else "gatk --java-options -Xmx" + memory + "G"

    command {
        set -e -o pipefail
        ~{preCommand}
        ~{toolCommand} \
        SplitNCigarReads \
        -I ~{inputBam} \
        -R ~{refFasta} \
        -O ~{outputBam} \
        -L ~{sep=' -L ' intervals}
    }

    output {
        File bam = outputBam
        File bamIndex = sub(outputBam, "\.bam$", ".bai")
    }

    runtime {
        memory: ceil(memory * memoryMultiplier)
    }
}
