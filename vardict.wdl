version 1.0

# Copyright (c) 2017 Leiden University Medical Center
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import "common.wdl"

task VarDict {
    input {
        String tumorSampleName
        File tumorBam
        File tumorBamIndex
        File referenceFasta
        File referenceFastaFai
        File bedFile
        String outputVcf
        Boolean outputCandidateSomaticOnly = true
        Boolean outputAllVariantsAtSamePosition = true
        Float mappingQuality = 20
        Int minimumTotalDepth = 8
        Int minimumVariantDepth = 4
        Float minimumAlleleFrequency = 0.02
        Int chromosomeColumn = 1
        Int startColumn = 2
        Int endColumn = 3
        Int geneColumn = 4

        String? normalSampleName
        File? normalBam
        File? normalBamIndex

        String javaXmx = "16G"
        Int threads = 1
        String memory = "18G"
        Int timeMinutes = 300
        String dockerImage = "quay.io/biocontainers/vardict-java:1.5.8--1"
    }

    command {
        set -e -o pipefail
        export JAVA_OPTS="-Xmx~{javaXmx} -XX:ParallelGCThreads=1"
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
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        tumorSampleName: {description: "The name of the tumor/case sample.", category: "required"}
        tumorBam: {description: "The tumor/case sample's BAM file.", category: "required"}
        tumorBamIndex: {description: "The index for the tumor/case sample's BAM file.", category: "required"}
        referenceFasta: {description: "The reference fasta file.", category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        bedFile: {description: "A bed file describing the regions to operate on. These regions must be below 1e6 bases in size.", category: "required"}
        outputVcf: {description: "The location to write the output VCF file to.", category: "required"}
        outputCandidateSomaticOnly: {description: "Equivalent to var2vcf_paired.pl or var2vcf_valid.pl's `-M` flag.", category: "advanced"}
        outputAllVariantsAtSamePosition: {description: "Equivalent to var2vcf_paired.pl or var2vcf_valid.pl's `-A` flag.", category: "advanced"}
        mappingQuality: {description: "Equivalent to var2vcf_paired.pl or var2vcf_valid.pl's `-Q` option.", category: "advanced"}
        minimumTotalDepth: {description: "Equivalent to var2vcf_paired.pl or var2vcf_valid.pl's `-d` option.", category: "advanced"}
        minimumVariantDepth: {description: "Equivalent to var2vcf_paired.pl or var2vcf_valid.pl's `-v` option.", category: "advanced"}
        minimumAlleleFrequency: {description: "Equivalent to var2vcf_paired.pl or var2vcf_valid.pl's `-f` option.", category: "advanced"}
        chromosomeColumn: {description: "Equivalent to vardict-java's `-c` option.", category: "advanced"}
        startColumn: {description: "Equivalent to vardict-java's `-S` option.", category: "advanced"}
        endColumn: {description: "Equivalent to vardict-java's `-E` option.", category: "advanced"}
        geneColumn: {description: "Equivalent to vardict-java's `-g` option.", category: "advanced"}
        normalSampleName: {description: "The name of the normal/control sample.", category: "common"}
        normalBam: {description: "The normal/control sample's BAM file.", category: "common"}
        normalBamIndex: {description: "The normal/control sample's BAM file.", category: "common"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.", category: "advanced"}
        threads: {description: "The number of threads to use.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}
    }
}
