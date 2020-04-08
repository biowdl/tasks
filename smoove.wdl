version 1.0

# MIT License
#
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

task CallSV {
    input {
        File bamFile
        File bamIndex
        File referenceFasta
        File referenceFastaFai
        String sample
        String outputDir = "./smoove"

        String memory = "15G"
        String dockerImage = "quay.io/biocontainers/smoove:0.2.5--0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputDir})"
        smoove call \
        --outdir ~{outputDir} \
        --name ~{sample} \
        --fasta ~{referenceFasta} \
        ~{bamFile}
    }

    output {
        File smooveVcf = outputDir + "/" + sample + "-smoove.vcf.gz"
    }

    runtime {
        memory: memory
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        bamFile: {description: "The bam file to process.", category: "required"}
        bamIndex: {description: "The index bam file.", category: "required"}
        referenceFasta: {description: "The reference fasta file also used for mapping.", category: "required"}
        referenceFastaFai: {description: "Fasta index (.fai) file of the reference.", category: "required" }
        outputDir: {description: "The location the output VCF file should be written.", category: "common"}
        sample: {description: "The name of the sample.", category: "required"}
        memory: {description: "The memory required to run the programs.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}
    }
}
