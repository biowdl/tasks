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
        Array[File]+ bamFile
        Array[File]+ bamIndex
        File referenceFasta
        File referenceFastaFai
        String outputPath = "./delly/delly.bcf"

        File? genotypeBcf
        File? genotypeBcfIndex

        String memory = "15GiB"
        Int timeMinutes = 600
        String dockerImage = "quay.io/biocontainers/delly:1.1.6--ha41ced6_0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPath})"
        delly call \
        -o ~{outputPath} \
        -g ~{referenceFasta} \
        ~{"-v " + genotypeBcf} \
        ~{sep=" " bamFile}
    }

    output {
        File dellyBcf = outputPath
        File dellyBcfIndex = outputPath + ".csi"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        bamFile: {description: "The bam files to process.", category: "required"}
        bamIndex: {description: "The indexes for the bam files.", category: "required"}
        referenceFasta: {description: "The reference fasta file also used for mapping.", category: "required"}
        referenceFastaFai: {description: "Fasta index (.fai) file of the reference.", category: "required" }
        outputPath: {description: "The location the output BCF file should be written.", category: "common"}
        genotypeBcf: {description: "A BCF with SVs to get genotyped in the samples.", category: "advanced"}
        genotypeBcfIndex: {description: "The index for the genotype BCF file.", category: "advanced"}
        memory: {description: "The memory required to run the programs.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        dellyBcf: {description: "File containing structural variants."}
    }
}


task SomaticFilter {
    input {
        File dellyBcf
        File dellyBcfIndex
        Array[String]+ normalSamples
        Array[String]+ tumorSamples
        String outputPath = "./delly/delly_filter.bcf"

        String memory = "15GiB"
        Int timeMinutes = 300
        String dockerImage = "quay.io/biocontainers/delly:1.1.6--ha41ced6_0"
    }

    command <<<
        set -e
        mkdir -p "$(dirname ~{outputPath})"
        for SAMPLE in ~{sep=" " normalSamples}; do echo -e "${SAMPLE}\tcontrol" >> samples.tsv; done
        for SAMPLE in ~{sep=" " tumorSamples}; do echo -e "${SAMPLE}\ttumor" >> samples.tsv; done

        delly filter \
        -f somatic \
        -o ~{outputPath} \
        -s samples.tsv \
        ~{dellyBcf}
    >>>

    output {
        File filterBcf = outputPath
        File filterBcfIndex = outputPath + ".csi"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        dellyBcf: {description: "The BCF file produced by delly.", category: "required"}
        dellyBcfIndex: {description: "The index for the delly BCF file.", category: "required"}
        normalSamples: {description: "The names for the normal samples as used in the delly BCF file.", category: "required"}
        tumorSamples: {description: "The names for the tumor samples as used in the delly BCF file.", category: "required"}
        outputPath: {description: "The location the output BCF file should be written.", category: "common"}
        memory: {description: "The memory required to run the programs.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}
    }
}