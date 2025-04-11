version 1.0

# Copyright (c) 2020 Leiden University Medical Center
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

task Normalize {
    input {
        File inputVCF
        File inputVCFIndex
        File referenceFasta
        File referenceFastaFai
        Boolean ignoreMaskedRef = false
        String outputPath = "./vt/normalized_decomposed.vcf.gz"
        String? filterExpression

        Int compressionLevel = 1

        String memory = "4GiB"
        Int timeMinutes = 10 + ceil(size(inputVCF, "GiB") * 240)
        String dockerImage = "quay.io/biocontainers/vt:0.57721--h2419454_12"
    }

    command {
        set -eo pipefail
        mkdir -p "$(dirname ~{outputPath})"
        vt view -h \
        ~{"-f '" + filterExpression}~{true="'" false="" defined(filterExpression)} \
        ~{inputVCF} \
        | vt normalize - \
        -r ~{referenceFasta} \
        ~{true="-m " false="" ignoreMaskedRef} \
        | vt decompose -s - \
        | vt view - \
        -c ~{compressionLevel} \
        -o ~{outputPath} 
        vt index ~{outputPath}
    }

    output {
        File outputVcf = outputPath
        File outputVcfIndex = outputPath + ".tbi"
    }

    runtime {
        cpu: 1
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        inputVCF: {description: "The VCF file to process.", category: "required"}
        inputVCFIndex: {description: "The index of the VCF file to be processed.", category: "required"}
        referenceFasta: {description: "The reference fasta file which was also used for mapping.", category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        ignoreMaskedRef: {description: "Warns but does not exit when REF is inconsistent with masked reference sequence for non SNPs.", category: "advanced"}
        outputPath: {description: "The location the output VCF file should be written.", category: "common"}
        filterExpression: {description: "See https://genome.sph.umich.edu/wiki/Vt#Filters for valid expressions.", category: "common"}
        compressionLevel: {description: "Compression level for the out vcf.gz file.", category: "advanced"}

        memory: {description: "The memory required to run the programs.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        outputVcf: {description: "Normalized and decomposed VCF file."}
        outputVcfIndex: {description: "Index for normalized and decomposed VCF file."}
    }
}
