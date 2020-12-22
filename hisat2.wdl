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

task Hisat2 {
    input {
        File inputR1
        File? inputR2
        Array[File]+ indexFiles
        String outputBam
        String sample
        String library
        String readgroup
        String platform = "illumina"
        Boolean downstreamTranscriptomeAssembly = true
        String summaryFilePath = basename(outputBam, ".bam") + ".summary.txt"
        Int sortMemoryPerThreadGb = 2
        Int compressionLevel = 1

        Int? sortThreads

        Int threads = 4
        Int? memoryGb
        Int timeMinutes = 1 + ceil(size([inputR1, inputR2], "G") * 180 / threads)
        # quay.io/biocontainers/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1
        # is a combination of hisat2 and samtools hisat2=2.2.0 & samtools=1.10.
        String dockerImage = "quay.io/biocontainers/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:2880dd9d8ad0a7b221d4eacda9a818e92983128d-0"
    }

    # Samtools sort may block the pipe while it is writing data to disk.
    # This can lead to cpu underutilization.
    # 1 thread if threads is 1. For 2-4 threads 2 sort threads. 3 sort threads for 5-8 threads.
    Int estimatedSortThreads = if threads == 1 then 1 else 1 + ceil(threads / 4.0)
    Int totalSortThreads = select_first([sortThreads, estimatedSortThreads])
    Int estimatedMemoryGb = 1 + ceil(size(indexFiles, "G") * 1.2) + sortMemoryPerThreadGb * totalSortThreads

    command {
        set -e -o pipefail
        mkdir -p "$(dirname ~{outputBam})"
        hisat2 \
        -p ~{threads} \
        -x ~{sub(indexFiles[0], "\.[0-9]\.ht2", "")} \
        ~{true="-1" false="-U" defined(inputR2)} ~{inputR1} \
        ~{"-2" + inputR2} \
        --rg-id ~{readgroup} \
        --rg 'SM:~{sample}' \
        --rg 'LB:~{library}' \
        --rg 'PL:~{platform}' \
        ~{true="--dta" false="" downstreamTranscriptomeAssembly} \
        --new-summary \
        --summary-file ~{summaryFilePath} \
        | samtools sort \
        ~{"-@ " + totalSortThreads} \
        -m ~{sortMemoryPerThreadGb}G \
        -l ~{compressionLevel} \
        - \
        -o ~{outputBam}
    }

    output {
        File bamFile = outputBam
        File summaryFile = summaryFilePath
    }

    runtime {
        cpu: threads
        memory: "~{select_first([memoryGb, estimatedMemoryGb])}G"
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        inputR1: {description: "The first-/single-end FastQ file.", category: "required"}
        inputR2: {description: "The second-end FastQ file.", category: "common"}
        indexFiles: {description: "The hisat2 index files.", category: "required"}
        outputBam: {description: "The location the output BAM file should be written to.", category: "required"}
        sample: {description: "The sample id.", category: "required"}
        library: {description: "The library id.", category: "required"}
        readgroup: {description: "The readgroup id.", category: "required"}
        platform: {description: "The platform used for sequencing.", category: "advanced"}
        downstreamTranscriptomeAssembly: {description: "Equivalent to hisat2's `--dta` flag.", category: "advanced"}
        summaryFilePath: {description: "Where the summary file should be written.", category: "advanced"}
        sortMemoryPerThreadGb: {description: "The amount of memory for each sorting thread in gigabytes.", category: "advanced"}
        compressionLevel: {description: "The compression level of the output BAM.", category: "advanced"}
        sortThreads: {description: "The number of threads to use for sorting.", category: "advanced"}
        threads: {description: "The number of threads to use.", category: "advanced"}
        memoryGb: {description: "The amount of memory this job will use in gigabytes.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        bamFile: {description: "Output BAM file."}
        summaryFile: {description: "Alignment summary file."}
    }
}
