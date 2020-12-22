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

import "bwa.wdl"

task Mateclever {
    input {
        File fiteredBam
        File indexedFiteredBam
        BwaIndex bwaIndex
        File predictions
        String outputPath = "./clever"
        Int cleverMaxDelLength = 100000
        Int maxLengthDiff= 30
        Int maxOffset = 150

        Int threads = 10
        String memory = "15G"
        Int timeMinutes = 600
        String dockerImage = "quay.io/biocontainers/clever-toolkit:2.4--py36hcfe0e84_6"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPath})"
        echo ~{outputPath} ~{fiteredBam} ~{predictions} none > predictions.list
        mateclever \
        -T ~{threads} \
        -k \
        -f \
        -M ~{cleverMaxDelLength} \
        -z ~{maxLengthDiff} \
        -o ~{maxOffset} \
        ~{bwaIndex.fastaFile} \
        predictions.list \
        ~{outputPath}
    }

    output {
        File matecleverVcf = outputPath + "/deletions.vcf"
    }

    runtime {
        cpu: threads
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        fiteredBam: {description: "The bam file where sequences less than 30bp were removed.", category: "required"}
        indexedFiteredBam: {description: "The index of the filtered bam file.", category: "required"}
        bwaIndex: {description: "The BWA index files.", category: "required"}
        predictions: {description: "The predicted deletions (VCF) from clever.", category: "required"}
        outputPath: {description: "The location the output VCF file should be written.", category: "common"}
        cleverMaxDelLength: {description: "The maximum deletion length to look for in Clever predictions.", category: "advanced"}
        maxLengthDiff: {description: "The maximum length difference between split-read and read-pair deletion to be considered identical.", category: "advanced"}
        maxOffset: {description: "The maximum center distance between split-read and read-pair deletion to be considered identical.", category: "advanced"}
        threads: {description: "The the number of threads required to run a program.", category: "advanced"}
        memory: {description: "The memory required to run the programs.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        matecleverVcf: {description: "VCF with additional mateclever results."}
    }
}

task Prediction {
    input {
        File bamFile
        File bamIndex
        BwaIndex bwaIndex
        String outputPath = "./clever"

        Int threads = 10
        String memory = "55G"
        Int timeMinutes = 480
        String dockerImage = "quay.io/biocontainers/clever-toolkit:2.4--py36hcfe0e84_6"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPath})"
        clever \
        -T ~{threads} \
        --use_mapq \
        --sorted \
        -f \
        ~{bamFile} \
        ~{bwaIndex.fastaFile} \
        ~{outputPath}
    }

    output {
        File predictions = outputPath + "/predictions.vcf"
    }

    runtime {
        cpu: threads
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        bamFile: {description: "The bam file to process.", category: "required"}
        bamIndex: {description: "The index bam file.", category: "required"}
        bwaIndex: {description: "The BWA index files.", category: "required"}
        outputPath: {description: "The location the output VCF file should be written.", category: "common"}
        threads: {description: "The the number of threads required to run a program.", category: "advanced"}
        memory: {description: "The memory required to run the programs.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        predictions: {description: "The predicted deletions (VCF) from clever."}
    }
}
