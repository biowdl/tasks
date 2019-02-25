version 1.0

import "common.wdl"

task VarDict {
    input {
        String? installDir
        Boolean useJavaVersion = true

        String tumorSampleName
        IndexedBamFile tumorBam
        String? normalSampleName
        IndexedBamFile? normalBam
        Reference reference
        File bedFile
        String outputVcf

        Int chromosomeColumn = 1
        Int startColumn = 2
        Int endColumn = 3
        Int geneColumn = 4

        String? preCommand
        Int memory = 8
        Float memoryMultiplier = 2.0
    }

    String normalArg = if (defined(normalBam))
        then "|" + select_first([normalBam]).file
        else ""

    String toolCommand = if defined(installDir)
        then installDir + "/VarDict"
        else if useJavaVersion
            then "vardict-java"
            else "vardict"

    command {
        set -e -o pipefail
        export JAVA_OPTS="-Xmx~{memory}G"
        ~{preCommand}
        ~{toolCommand} \
        -G ~{reference.fasta} \
        -N ~{tumorSampleName} \
        -b "~{tumorBam.file}~{normalArg}" \
        ~{true="" false="-z" defined(normalBam)} \
        -c ~{chromosomeColumn} \
        -S ~{startColumn} \
        -E ~{endColumn} \
        -g ~{geneColumn} \
        ~{bedFile} | \
        ~{installDir + "/"}~{true="testsomatic.R" false="teststrandbias.R" defined(normalBam)} | \
        ~{installDir + "/"}~{true="var2vcf_paired.pl"
            false="var2vcf_valid.pl" defined(normalBam)} \
        -N "~{tumorSampleName}~{"|" + normalSampleName}" \
        ~{true="" false="-E" defined(normalBam)} | \
        bgzip -c > ~{outputVcf}
        tabix -p vcf ~{outputVcf}
    }

    output {
        IndexedVcfFile vcfFile = object {
          file: outputVcf,
          index: outputVcf + ".tbi"
        }
    }

    runtime {
        memory: ceil(memory * memoryMultiplier)
    }
}
