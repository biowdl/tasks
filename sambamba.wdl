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

task Flagstat {
    input {
        File inputBam
        File inputBamIndex
        String outputPath = "./flagstat.txt"

        Int threads = 2
        String memory = "8G"
        Int timeMinutes = 320
        String dockerImage = "quay.io/biocontainers/sambamba:0.7.1--h148d290_2"
    }

    command {
        sambamba flagstat \
        -t ~{threads} \
        ~{inputBam} \
        > ~{outputPath}
    }

    output {
        File stats = outputPath
    }

    runtime {
        cpu: threads
        memory: memory
        time_minutes: timeMinutes # !UnknownRuntimeKey
        docker: dockerImage
    }

    parameter_meta {
        inputBam: {description: "The input BAM file.", category: "required"}
        inputBamIndex: {description: "The index for the BAM file.", category: "required"}
        outputPath: {description: "The path to write the ouput to.", category: "required"}

        threads: {description: "The number of threads that will be used for this task.", category: "advanced"}
        memory: {description: "The amount of memory available to the job.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}
    }
}


task Markdup {
    input {
        Array[File] inputBams
        String outputPath
        Int compressionLevel = 1
        # sortBufferSize and ioBufferSize taken from markdup defaults as of sambamba 0.7.1.
        Int sortBufferSize = 4096
        Int ioBufferSize = 128
        Boolean removeDuplicates = false

        Int? hashTableSize
        Int? overFlowListSize

        # Sambamba scales like this: 1 thread is fully utilized (1).
        # 2 threads 1.8 utilized. 3 -> 2.4, 4-> 2.7.
        # 2 threads reduces wall clock time by more than 40%.
        Int threads = 2
        # According to the manual sambamba markdup uses the sortbufferSize + 2 times the ioBuffer size.
        # Added 8192 mb as a margin of safety. Real life use with this setting uses 2.7 GiB.
        Int memoryMb = 8192 + sortBufferSize + 2 * ioBufferSize
        # Time minute calculation does not work well for higher number of threads.
        Int timeMinutes = 1 + ceil(size(inputBams, "G") * 25) / threads
        String dockerImage = "quay.io/biocontainers/sambamba:0.7.1--h148d290_2"
    }

    String bamIndexPath = sub(outputPath, "\.bam$", ".bai")

    command {
        set -e 
        mkdir -p "$(dirname ~{outputPath})"
        sambamba markdup \
        --nthreads ~{threads} \
        -l ~{compressionLevel} \
        ~{true="-r" false="" removeDuplicates} \
        ~{"--hash-table-size " + hashTableSize} \
        ~{"--overflow-list-size " + overFlowListSize} \
        ~{"--sort-buffer-size " + sortBufferSize} \
        ~{"--io-buffer-size " + ioBufferSize} \
        ~{sep=' ' inputBams} ~{outputPath}
        # sambamba creates an index for us.
        mv ~{outputPath}.bai ~{bamIndexPath}
    }

    output {
        File outputBam = outputPath
        File outputBamIndex = bamIndexPath
    }

    runtime {
        cpu: threads
        memory: "~{memoryMb}M"
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        inputBams: {description: "The input BAM files.", category: "required"}
        outputPath: {description: "Output directory path + output file.", category: "required"}
        compressionLevel: {description: "Compression level from 0 (uncompressed) to 9 (best).", category: "advanced"}
        sortBufferSize: {description: "The amount of mb allocated to the sort buffer.", category: "advanced"}
        ioBufferSize: {description: "The amount of mb allocated to each IO buffer. Sambamba uses two IO buffers.", category: "advanced"}
        removeDuplicates: {description: "Whether to remove the duplicates (instead of only marking them).", category: "advanced"}
        hashTableSize: {description: "Sets sambamba's hash table size.", category: "advanced"}
        overFlowListSize: {description: "Sets sambamba's overflow list size.", category: "advanced"}
        threads: {description: "The number of threads that will be used for this task.", category: "advanced"}
        memoryMb: {description: "The amount of memory available to the job in megabytes.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        outputBam: {description: "Sorted BAM file."}
        outputBamIndex: {description: "Sorted BAM file index."}
    }
}

task Sort {
    input {
        File inputBam
        String outputPath = basename(inputBam, "\.bam") + ".sorted.bam"
        Boolean sortByName = false
        Int compressionLevel = 1

        Int memoryPerThreadGb = 4
        Int threads = 1
        Int memoryGb = 1 + threads * memoryPerThreadGb
        Int timeMinutes = 1 + ceil(size(inputBam, "G") * 3)
        String dockerImage = "quay.io/biocontainers/sambamba:0.7.1--h148d290_2"
    }

    # Select first needed as outputPath is optional input (bug in cromwell).
    String bamIndexPath = sub(select_first([outputPath]), "\.bam$", ".bai")

    command {
        set -e
        mkdir -p "$(dirname ~{outputPath})"
        sambamba sort \
        -l ~{compressionLevel} \
        ~{true="-n" false="" sortByName} \
        ~{"--nthreads " + threads} \
        -m ~{memoryPerThreadGb}G \
        -o ~{outputPath} \
        ~{inputBam}
        # sambamba creates an index for us.
        mv ~{outputPath}.bai ~{bamIndexPath}
    }

    output {
        File outputBam = outputPath
        File outputBamIndex = bamIndexPath
    }

    runtime {
        cpu: threads
        memory: "~{memoryGb}G"
        docker: dockerImage
        time_minutes: timeMinutes
    }

    parameter_meta {
        # inputs
        inputBam: {description: "The input SAM file.", category: "required"}
        outputPath: {description: "Output directory path + output file.", category: "required"}
        sortByName: {description: "Sort the inputBam by read name instead of position.", category: "advanced"}
        compressionLevel: {description: "Compression level from 0 (uncompressed) to 9 (best).", category: "advanced"}
        memoryPerThreadGb: {description: "The amount of memory used per sort thread in gigabytes.", category: "advanced"}
        threads: {description: "The number of threads that will be used for this task.", category: "advanced"}
        memoryGb: {description: "The amount of memory available to the job in gigabytes.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        outputBam: {description: "Sorted BAM file."}
        outputBamIndex: {description: "Sorted BAM file index."}
    }
}
