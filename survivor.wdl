version 1.0

# MIT License
#
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

task Merge {
    input{
        Array[File] filePaths
        Int breakpointDistance = 1000
        Int suppVecs = 2
        Boolean svType = true
        Boolean strandType = true
        Boolean distanceBySvSize = false
        Int minSize = 30
        String outputPath = "./survivor/merged.vcf"
        String memory = "24G"
        Int timeMinutes = 60
        String dockerImage = "quay.io/biocontainers/survivor:1.0.6--h6bb024c_0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPath})"
        echo '~{sep="\n" filePaths}' > fileList
        SURVIVOR merge \
        fileList \
        ~{breakpointDistance} \
        ~{suppVecs} \
        ~{true=1 false=0 svType} \
        ~{true=1 false=0 strandType} \
        ~{true=1 false=0 distanceBySvSize} \
        ~{minSize} \
        ~{outputPath}
    }

    output {
        File mergedVcf = outputPath
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        filePaths: {description: "An array of VCF files (predictions) to be merged by SURVIVOR", category: "required"}
        breakpointDistance: {description: "The distance between pairwise breakpoints between SVs", category: "advanced"}
        suppVecs: {description: "The minimum number of SV callers to support the merging", category: "advanced"}
        svType: {description: "A boolean to include the type SV to be merged", category: "advanced"}
        strandType: {description: "A boolean to include strand type of an SV to be merged", category: "advanced"}
        distanceBySvSize: {description: "A boolean to predict the pairwise distance between the SVs based on their size", category: "advanced"}
        minSize: {description: "The mimimum size of SV to be merged", category: "advanced"}
        outputPath: {description: "The location the output VCF file should be written.", category: "common"}
        memory: {description: "The memory required to run the programs", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}
    }
}
