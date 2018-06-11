# Apply Base Quality Score Recalibration (BQSR) model
task ApplyBQSR {
    String? preCommand
    File gatkJar
    File inputBam
    File inputBamIndex
    String outputBamPath
    File recalibrationReport
    Array[File]+ sequenceGroupInterval
    File refDict
    File refFasta
    File refFastaIndex
    Int? compressionLevel

    Float? memory
    Float? memoryMultiplier

    Int mem = ceil(select_first([memory, 4.0]))
    command {
        set -e -o pipefail
        ${preCommand}
        java ${"-Dsamjdk.compression_level=" + compressionLevel} \
        -Xms${mem}G -jar ${gatkJar} \
          ApplyBQSR \
          --create-output-bam-md5 \
          --add-output-sam-program-record \
          -R ${refFasta} \
          -I ${inputBam} \
          --use-original-qualities \
          -O ${outputBamPath} \
          -bqsr ${recalibrationReport} \
          --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
          -L ${sep=" -L " sequenceGroupInterval}
    }

    output {
        File recalibrated_bam = outputBamPath
        File recalibrated_bam_checksum = outputBamPath + ".md5"
    }

    runtime {
        memory: ceil(mem * select_first([memoryMultiplier, 3.0]))
    }
}

# Generate Base Quality Score Recalibration (BQSR) model
task BaseRecalibrator {
    String? preCommand
    File gatkJar
    File inputBam
    File inputBamIndex
    String recalibrationReportPath
    Array[File]+ sequenceGroupInterval
    Array[File]+ knownIndelsSitesVCFs
    Array[File]+ knownIndelsSitesIndices
    File refDict
    File refFasta
    File refFastaIndex

    Float? memory
    Float? memoryMultiplier

    Int mem = ceil(select_first([memory, 4.0]))
    command {
        set -e -o pipefail
        ${preCommand}
        java -Xms${mem}G -jar ${gatkJar} \
          BaseRecalibrator \
          -R ${refFasta} \
          -I ${inputBam} \
          --use-original-qualities \
          -O ${recalibrationReportPath} \
          --known-sites ${sep=" --known-sites " knownIndelsSitesVCFs} \
          -L ${sep=" -L " sequenceGroupInterval}
    }

    output {
        File recalibrationReport = recalibrationReportPath
    }

    runtime {
        memory: ceil(mem * select_first([memoryMultiplier, 3.0]))
    }
}

task CombineGVCFs {
    String? preCommand
    Array[File]+ gvcfFiles
    Array[File]+ gvcfFileIndexes
    Array[File]+ intervals

    String outputPath

    String gatkJar

    File refFasta
    File refFastaIndex
    File refDict

    Int? compressionLevel
    Float? memory
    Float? memoryMultiplier

    Int mem = ceil(select_first([memory, 4.0]))
    command {
        set -e -o pipefail
        ${preCommand}

        if [ ${length(gvcfFiles)} -gt 1 ]; then
            java ${"-Dsamjdk.compression_level=" + compressionLevel} \
            -Xmx${mem}G -jar ${gatkJar} \
             CombineGVCFs \
             -R ${refFasta} \
             -O ${outputPath} \
             -V ${sep=' -V ' gvcfFiles} \
             -L ${sep=' -L ' intervals}
        else # TODO this should be handeled in wdl
            ln -sf ${select_first(gvcfFiles)} ${outputPath}
            ln -sf ${select_first(gvcfFileIndexes)} ${outputPath}.tbi
        fi
    }

    output {
        File outputGVCF = outputPath
        File outputGVCFindex = outputPath + ".tbi"
    }

    runtime {
        memory: ceil(mem * select_first([memoryMultiplier, 3.0]))
    }
}

# Combine multiple recalibration tables from scattered BaseRecalibrator runs
task GatherBqsrReports {
    String? preCommand
    String gatkJar
    Array[File] inputBQSRreports
    String outputReportPath

    Float? memory
    Float? memoryMultiplier

    Int mem = ceil(select_first([memory, 4.0]))
    command {
        set -e -o pipefail
        ${preCommand}
        java -Xms${mem}G -jar ${gatkJar} \
        GatherBQSRReports \
        -I ${sep=' -I ' inputBQSRreports} \
        -O ${outputReportPath}
    }

    output {
        File outputBQSRreport = outputReportPath
    }

    runtime {
        memory: ceil(mem * select_first([memoryMultiplier, 3.0]))
    }
}

task GenotypeGVCFs {
    String? preCommand
    File gvcfFiles
    File gvcfFileIndexes
    Array[File]+ intervals

    String outputPath

    String gatkJar

    File refFasta
    File refFastaIndex
    File refDict

    File dbsnpVCF
    File dbsnpVCFindex

    Int? compressionLevel
    Float? memory
    Float? memoryMultiplier

    Int mem = ceil(select_first([memory, 4.0]))
    command {
        set -e -o pipefail
        ${preCommand}

        java ${"-Dsamjdk.compression_level=" + compressionLevel} \
        -Xmx${mem}G -jar ${gatkJar} \
         GenotypeGVCFs \
         -R ${refFasta} \
         -O ${outputPath} \
         -D ${dbsnpVCF} \
         -G StandardAnnotation \
         --only-output-calls-starting-in-intervals \
         -new-qual \
         -V ${gvcfFiles} \
         -L ${sep=' -L ' intervals}
    }

    output {
        File outputVCF = outputPath
        File outputVCFindex = outputPath + ".tbi"
    }

    runtime{
        memory: ceil(mem * select_first([memoryMultiplier, 3.0]))
    }
}

# Call variants on a single sample with HaplotypeCaller to produce a GVCF
task HaplotypeCallerGvcf {
    String? preCommand
    Array[File]+ inputBams
    Array[File]+ inputBamsIndex
    Array[File]+ intervalList
    String gvcfPath
    File refDict
    File refFasta
    File refFastaIndex
    Float? contamination
    Int? compressionLevel
    String gatkJar

    Float? memory
    Float? memoryMultiplier

    Int mem = ceil(select_first([memory, 4.0]))
    command {
        set -e -o pipefail
        ${preCommand}
        java ${"-Dsamjdk.compression_level=" + compressionLevel} \
        -Xmx${mem}G -jar ${gatkJar} \
          HaplotypeCaller \
          -R ${refFasta} \
          -O ${gvcfPath} \
          -I ${sep=" -I " inputBams} \
          -L ${sep=' -L ' intervalList} \
          -contamination ${default=0 contamination} \
          -ERC GVCF
    }

    output {
        File outputGVCF = gvcfPath
        File outputGVCFindex = gvcfPath + ".tbi"
    }

    runtime {
        memory: ceil(mem * select_first([memoryMultiplier, 3.0]))
    }
}

task SplitNCigarReads {
    String? preCommand

    File inputBam
    File inputBamIndex
    File refFasta
    File refFastaIndex
    File refDict
    String outputBam
    String gatkJar
    Array[File]+ intervals

    Float? memory
    Float? memoryMultiplier

    Int mem = ceil(select_first([memory, 4.0]))
    command {
        set -e -o pipefail
        ${preCommand}
        java -Xms${mem}G -jar ${gatkJar} \
        SplitNCigarReads \
        -I ${inputBam} \
        -R ${refFasta} \
        -O ${outputBam} \
        -L ${sep=' -L ' intervals}
    }

    output {
        File bam = outputBam
        File bam_index = sub(outputBam, "\\.bam$", ".bai")
    }

    runtime {
        memory: ceil(mem * select_first([memoryMultiplier, 3.0]))
    }
}
