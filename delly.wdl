version 1.0

# Copyright (c) 2018 Leiden University Medical Center
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
        String outputPath = "./delly/delly.bcf"

        String memory = "15GiB"
        Int timeMinutes = 300
        String dockerImage = "quay.io/biocontainers/delly:0.8.1--h4037b6b_1"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPath})"
        delly call \
        -o ~{outputPath} \
        -g ~{referenceFasta} \
        ~{bamFile}
    }

    output {
        File dellyBcf = outputPath
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        bamFile: {description: "The bam file to process.", category: "required"}
        bamIndex: {description: "The index bam file.", category: "required"}
        referenceFasta: {description: "The reference fasta file also used for mapping.", category: "required"}
        referenceFastaFai: {description: "Fasta index (.fai) file of the reference.", category: "required" }
        outputPath: {description: "The location the output VCF file should be written.", category: "common"}
        memory: {description: "The memory required to run the programs.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        dellyBcf: {description: "File containing structural variants."}
    }
}
