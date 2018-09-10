version 1.0

import "common.wdl"

# Apply Base Quality Score Recalibration (BQSR) model
task ApplyBQSR {
    input {
        String? preCommand
        File? gatkJar
        IndexedBamFile inputBam
        String outputBamPath
        File recalibrationReport
        Array[File]+ sequenceGroupInterval
        Reference reference

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
        IndexedBamFile recalibrated_bam = {
            "file": outputBamPath,
            "index": sub(outputBamPath, "\.bam$", ".bai")
        }

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
        IndexedBamFile inputBam
        String recalibrationReportPath
        Array[File]+ sequenceGroupInterval
        Array[IndexedVcfFile]? knownIndelsSitesVCFs
        IndexedVcfFile? dbsnpVCF
        Reference reference
        Int memory = 4
        Float memoryMultiplier = 3.0
    }

    Array[File]+ knownIndelsSitesVCFsArg = flatten([
        select_first([knownIndelsSitesVCFs.file, []]),
        select_all([dbsnpVCF.file])
    ])

    String toolCommand = if defined(gatkJar)
        then "java -Xmx" + memory + "G -jar " + gatkJar
        else "gatk --java-options -Xmx" + memory + "G"

    command {
        set -e -o pipefail
        ~{preCommand}
        ~{toolCommand} \
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
        memory: ceil(memory * memoryMultiplier)
    }
}

task CombineGVCFs {
    input {
        String? preCommand
        Array[IndexedVcfFile]+ gvcfFiles
        Array[File]+ intervals

        String outputPath

        String? gatkJar

        Reference reference

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
             -R ~{reference.fasta} \
             -O ~{outputPath} \
             -V ~{sep=' -V ' gvcfFiles} \
             -L ~{sep=' -L ' intervals}
        else # TODO this should be handeled in wdl
            ln -sf ~{gvcfFiles[0].file} ~{outputPath}
            ln -sf ~{gvcfFiles[0].index} ~{outputPath}.tbi
        fi
    }

    output {
        IndexedVcfFile outputVCF = {
            "file": outputPath,
            "index": outputPath + ".tbi"
        }
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
        Array[IndexedVcfFile] gvcfFiles
        Array[File]+ intervals

        String outputPath

        String? gatkJar

        Reference reference

        IndexedVcfFile? dbsnpVCF

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
        -R ~{reference.fasta} \
        -O ~{outputPath} \
        ~{"-D " + dbsnpVCF.file} \
        -G StandardAnnotation \
        --only-output-calls-starting-in-intervals \
        -new-qual \
        -V ~{sep=' -V ' gvcfFiles.file} \
        -L ~{sep=' -L ' intervals}
    }

    output {
        IndexedVcfFile outputVCF = {
            "file": outputPath,
            "index": outputPath + ".tbi"
        }
    }

    runtime{
        memory: ceil(memory * memoryMultiplier)
    }
}

# Call variants on a single sample with HaplotypeCaller to produce a GVCF
task HaplotypeCallerGvcf {
    input {
        String? preCommand
        Array[IndexedBamFile]+ inputBams
        Array[File]+ intervalList
        String gvcfPath
        Reference reference
        Float contamination = 0.0
        String? gatkJar

        IndexedVcfFile? dbsnpVCF

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
        -R ~{reference.fasta} \
        -O ~{gvcfPath} \
        -I ~{sep=" -I " inputBams.file} \
        -L ~{sep=' -L ' intervalList} \
        ~{"-D " + dbsnpVCF.file} \
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
        memory: ceil(memory * memoryMultiplier)
    }
}

task MuTect2 {
    input {
        String? preCommand

        Array[IndexedBamFile]+ inputBams
        Reference reference
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
        -R ~{reference.fasta} \
        -I ~{sep=" -I " inputBams.file} \
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
        memory: ceil(memory * memoryMultiplier)
    }
}

task SplitNCigarReads {
    input {
        String? preCommand

        IndexedBamFile inputBam
        Reference reference
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
        memory: ceil(memory * memoryMultiplier)
    }
}
