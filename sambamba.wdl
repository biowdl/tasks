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

task Sort {
    input {
        File inputBam
        String outputPath = basename(inputBam, "\.bam") + ".sorted.bam"
        Boolean sortByName = false
        Int compressionLevel = 1
        Int threads = 1
        Int memoryPerThreadGb = 4
        Int memoryGb = 1 + threads * memoryPerThreadGb
        String dockerImage = "quay.io/biocontainers/sambamba:0.7.1--h148d290_2"
        Int timeMinutes = 1 + ceil(size(inputBam, "G") * 2)
    }

    # Select first needed as outputPath is optional input. (bug in cromwell)
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
        sambamba index \
        ~{"--nthreads " + threads} \
         ~{outputPath} ~{bamIndexPath}
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
        memoryGb: {description: "The amount of memory available to the job in gigabytes.", category: "advanced"}
        memoryPerThreadGb: {description: "The amount of memory used per sort thread in gigabytes", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}
        threads: {description: "The number of additional threads that will be used for this task.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        # outputs
        outputBam: {description: "Sorted BAM file."}
    }
}