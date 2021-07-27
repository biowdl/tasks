version 1.0

# Copyright (c) 2021 Leiden University Medical Center
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

task Peach {
    input {
        File transcriptTsv
        File germlineVcf
        File germlineVcfIndex
        File tumorName
        File normalName
        String outputDir = "./peach"
        File panelJson

        String memory = "8G"
        String dockerImage = "quay.io/biowdl/peach:v1.0"
        Int timeMinutes = 20
    }

    command {
        peach \
        --recreate_bed \
        --transcript_tsv ~{transcriptTsv} \
        ~{germlineVcf} \
        ~{tumorName} \
        ~{normalName} \
        1.0 \
        ~{outputDir} \
        ~{panelJson} \
        vcftools
    }

    output {
        File callsTsv = "~{outputDir}/~{tumorName}.peach.calls.tsv"
        File filteredVcf = "~{outputDir}/~{tumorName}.peach.filtered.vcf"
        File genotypeTsv = "~{outputDir}/~{tumorName}.peach.genotype.tsv"
        Array[File] peachFiles = [callsTsv, filterVcf, genotypeTsv]
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes # !UnknownRuntimeKey
        docker: dockerImage
    }

    parameter_meta {
        transcriptTsv: {description: "A tsv file describing transcripts.", category: "required"}
        germlineVcf: {description: "The germline VCF file from hmftools' purple.", category: "required"}
        germlineVcfIndex: {description: "The germline VCF's index.", category: "required"}
        tumorName: {description: "The name of the tumor sample.", category: "required"}
        normalName: {description: "The name of the normal sample", category: "required"}
        outputDir: {description: "The directory the ouput should be written to.", category: "required"}
        panelJson: {description: "A JSON describing the panel.", category: "required"}

        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        memory: {description: "The amount of memory available to the job.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}
    }
}