version 1.0

import "common.wdl"

task VarDict {
    input {
        String tumorSampleName
        File tumorBam
        File tumorBamIndex
        String? normalSampleName
        File? normalBam
        File? normalBamIndex
        File referenceFasta
        File referenceFastaFai
        File bedFile
        String outputVcf

        Int chromosomeColumn = 1
        Int startColumn = 2
        Int endColumn = 3
        Int geneColumn = 4

        Boolean outputCandidateSomaticOnly = true
        Boolean outputAllVariantsAtSamePosition = true
        Float mappingQuality = 20
        Int minimumTotalDepth = 8
        Int minimumVariantDepth = 4
        Float minimumAlleleFrequency = 0.02

        Int threads = 1
        String memory = "40G"
        String javaXmx = "16G"
        String dockerImage = "quay.io/biocontainers/vardict-java:1.5.8--1"
    }

    command {
        set -e -o pipefail
        export JAVA_OPTS="-Xmx~{javaXmx}"
        vardict-java \
        ~{"-th " + threads} \
        -G ~{referenceFasta} \
        -N ~{tumorSampleName} \
        -b "~{tumorBam}~{"|" + normalBam}" \
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
        ~{true="-M" false="" outputCandidateSomaticOnly} \
        ~{true="-A" false="" outputAllVariantsAtSamePosition} \
        -Q ~{mappingQuality} \
        -d ~{minimumTotalDepth} \
        -v ~{minimumVariantDepth} \
        -f ~{minimumAlleleFrequency} \
        > ~{outputVcf}
    }

    output {
        File vcfFile = outputVcf
    }

    runtime {
        cpu: threads + 2
        memory: memory
        docker: dockerImage
    }
}
