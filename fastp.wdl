version 1.0

# MIT License
#
# Copyright (c) 2022 Leiden University Medical Center
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

task Fastp {
    input {
        File read1
        File read2
        String outputPathR1
        String outputPathR2
        String htmlPath
        String jsonPath

        Int compressionLevel = 1
        Boolean correction = false
        Int lengthRequired = 15
        Int? split
        Boolean performAdapterTrimming = true
        
        Int threads = 4
        String memory = "20GiB"
        Int timeMinutes = 1 + ceil(size([read1, read2], "G")  * 7.0 / threads)
        String dockerImage = "quay.io/biocontainers/fastp:0.23.2--h5f740d0_3"
    }

    String outputDirR1 = sub(outputPathR1, basename(outputPathR1), "")
    String outputDirR2 = sub(outputPathR2, basename(outputPathR2), "")

    command <<<
        set -e 
        mkdir -p $(dirname ~{outputPathR1})
        mkdir -p $(dirname ~{outputPathR2})
        mkdir -p $(dirname ~{htmlPath})
        mkdir -p $(dirname ~{jsonPath})

        # predict output paths
        seq 1 ~{if defined(split) then split else "2"} | awk '{print "~{outputDirR1}/"$0".~{basename(outputPathR1)}"}' > r1_paths
        seq 1 ~{if defined(split) then split else "2"} | awk '{print "~{outputDirR2}/"$0".~{basename(outputPathR2)}"}' > r2_paths
        fastp \
        -i ~{read1} \
        ~{"-I " + read2} \
        -o ~{outputPathR1} \
        ~{"-O " + outputPathR2} \
        -h ~{htmlPath} \
        -j ~{jsonPath} \
        -z ~{compressionLevel} \
        ~{if correction then "--correction" else ""} \
        --length_required ~{lengthRequired} \
        --thread ~{threads} \
        ~{"--split " + split} \
        ~{if defined(split) then "-d 0" else ""} \
        ~{if performAdapterTrimming then "" else "--disable_adapter_trimming"}
    >>>

    output {
        File htmlReport = htmlPath
        File jsonReport = jsonPath
        Array[File] clippedR1 = if defined(split) then read_lines("r1_paths") else [outputPathR1]
        Array[File] clippedR2 = if defined(split) then read_lines("r2_paths") else [outputPathR2]
    }

    runtime {
        cpu: threads
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        read1: {description: "The R1 fastq file.", category: "required"}
        read2: {description: "The R2 fastq file.", category: "required"}
        outputPathR1: {description: "The output path for the R1 file.", category: "required"}
        outputPathR2: {description: "The output path for the R2 file.", category: "required"}
        htmlPath: {description: "The path to write the html report to.", category: "required"}
        jsonPath: {description: "The path to write the json report to.", category: "required"}
        compressionLevel: {description: "The compression level to use for the output.", category: "advanced"}
        correction: {description: "Whether or not to apply overlap based correction.", category: "advanced"}
        lengthRequired: {description: "The minimum read length.", category: "advanced"}
        split: {description: "The number of chunks to split the files into.", category: "common"}
        performAdapterTrimming: {description: "Whether adapter trimming should be performed or not.", category: "advanced"}
        threads: {description: "The number of threads to use.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}
    }
}