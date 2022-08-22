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

task ExtractSigPredictHRD {
    input {
        String outputDir = "."
        String sampleName
        File snvIndelVcf
        File snvIndelVcfIndex
        File svVcf
        File svVcfIndex
        Boolean hg38 = false

        String memory = "3GiB"
        Int timeMinutes = 10
        String dockerImage = "quay.io/biowdl/chord-mutsigextractor:2.00_1.14"
    }

    command {
        extractSigPredictHRD.R \
        ~{outputDir} \
        ~{sampleName} \
        ~{snvIndelVcf} \
        ~{svVcf} \
        ~{if hg38 then "RG_38" else "RG_37"}
    }

    output {
        File chordPrediction = "~{outputDir}/~{sampleName}_chord_prediction.txt"
        File chordSignatures = "~{outputDir}/~{sampleName}_chord_signatures.txt"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes # !UnknownRuntimeKey
        docker: dockerImage
    }

    parameter_meta {
        outputDir: {description: "The directory the outout will be written to.", category: "required"}
        sampleName: {description: "The name of the sample.", category: "required"}
        snvIndelVcf: {description: "A VCF file with SNVs and indels.", category: "required"}
        snvIndelVcfIndex: {description: "The index for the SNV/indel VCF file.", category: "required"}
        svVcf: {description: "A VCF file with SVs.", category: "required"}
        svVcfIndex: {description: "The index for the SV VCF file.", category: "required"}

        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}