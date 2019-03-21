version 1.0

import "common.wdl"

task VarDict {
    input {
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

        Int threads = 1
        Int memory = 16
        Float memoryMultiplier = 2.5
        String dockerTag = "1.5.8--1"

        File? doNotDefineThis #FIXME
    }

    File? normalBamFile = if defined(normalBam)
        then select_first([normalBam]).file
        else doNotDefineThis

    command {
        set -e -o pipefail
        export JAVA_OPTS="-Xmx~{memory}G"
        vardict-java \
        ~{"-th " + threads} \
        -G ~{reference.fasta} \
        -N ~{tumorSampleName} \
        -b "~{tumorBam.file}~{"|" + normalBamFile}" \
        ~{true="" false="-z" defined(normalBam)} \
        -c ~{chromosomeColumn} \
        -S ~{startColumn} \
        -E ~{endColumn} \
        -g ~{geneColumn} \
        ~{bedFile} | \
        ~{true="testsomatic.R" false="teststrandbias.R" defined(normalBam)} | \
        ~{true="var2vcf_paired.pl" false="var2vcf_valid.pl" defined(normalBam)} \
        -N "~{tumorSampleName}~{"|" + normalSampleName}" \
        ~{true="" false="-E" defined(normalBam)} \
        > ~{outputVcf}
    }

    output {
        File vcfFile = outputVcf
    }

    runtime {
        cpu: threads + 2
        memory: ceil(memory * memoryMultiplier)
        docker: "quay.io/biocontainers/vardict-java:" + dockerTag
    }
}
