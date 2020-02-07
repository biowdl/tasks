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

    parameter_meta {
        tumorSampleName: {description: "The name of the tumor/case sample.", category: "required"}
        tumorBam: {description: "The tumor/case sample's BAM file.", category: "required"}
        tumorBamIndex: {description: "The index for the tumor/case sample's BAM file.", category: "required"}
        normalSampleName: {description: "The name of the normal/control sample.", category: "common"}
        normalBam: {description: "The normal/control sample's BAM file.", category: "common"}
        normalBamIndex: {description: "The normal/control sample's BAM file.", category: "common"}
        referenceFasta: {description: "The reference fasta file.", category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        bedFile: {description: "A bed file describing the regions to operate on. These regions must be below 1e6 bases in size.", category: "required"}
        outputVcf: {description: "The location to write the output VCF file to.", category: "required"}
        chromosomeColumn: {description: "Equivalent to vardict-java's `-c` option.", category: "advanced"}
        startColumn: {description: "Equivalent to vardict-java's `-S` option.", category: "advanced"}
        endColumn: {description: "Equivalent to vardict-java's `-E` option.", category: "advanced"}
        geneColumn: {description: "Equivalent to vardict-java's `-g` option.", category: "advanced"}
        outputCandidateSomaticOnly: {description: "Equivalent to var2vcf_paired.pl or var2vcf_valid.pl's `-M` flag.", category: "advanced"}
        outputAllVariantsAtSamePosition: {description: "Equivalent to var2vcf_paired.pl or var2vcf_valid.pl's `-A` flag.", category: "advanced"}
        mappingQuality: {description: "Equivalent to var2vcf_paired.pl or var2vcf_valid.pl's `-Q` option.", category: "advanced"}
        minimumTotalDepth: {description: "Equivalent to var2vcf_paired.pl or var2vcf_valid.pl's `-d` option.", category: "advanced"}
        minimumVariantDepth: {description: "Equivalent to var2vcf_paired.pl or var2vcf_valid.pl's `-v` option.", category: "advanced"}
        minimumAlleleFrequency: {description: "Equivalent to var2vcf_paired.pl or var2vcf_valid.pl's `-f` option.", category: "advanced"}

        threads: {description: "The number of threads to use.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.",
                  category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}
