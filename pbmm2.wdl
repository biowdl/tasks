version 1.0

# Copyright (c) 2020 Leiden University Medical Center
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

task Mapping {
    input {
        String presetOption
        Boolean unmapped = false
        String sample
        String outputPrefix = "./~{sample}.align"
        File referenceMMI
        File queryFile

        Int sortMemoryGb = 1
        Int sortThreads = 2
        Int cores = 8
        String memory = "24GiB"
        # Slightly higher than minimap2 as compression level can not be set.
        Int timeMinutes = 1 + ceil(size(queryFile, "G") * 400 / cores)
        String dockerImage = "quay.io/biocontainers/pbmm2:1.17.0--h9ee0642_0"
    }

    # Use cores+sortThreads to set the number of threads. Internally pbmm2
    # allocates cores - sortThreads to alignment. This leads to underutilization
    # of the requested resources. Sorting uses very little CPU until the point
    # comes that the memory is full and the temporary file needs to be written.
    # At this point the alignment halts because the pipe is full.
    command {
        set -e 
        mkdir -p $(dirname ~{outputPrefix})
        pbmm2 align \
        --preset ~{presetOption} \
        --sort \
        ~{true="--unmapped" false="" unmapped} \
        --num-threads ~{cores + sortThreads} \
        --sort-memory ~{sortMemoryGb}G \
        --sort-threads ~{sortThreads} \
        ~{referenceMMI} \
        ~{queryFile} \
        --sample ~{sample} \
        ~{outputPrefix}.bam
    }

    output {
        File outputAlignmentFile = outputPrefix + ".bam"
        File outputIndexFile = outputPrefix + ".bam.bai"
    }

    runtime {
        cpu: cores
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        presetOption: {description: "This option applies multiple options at the same time.", category: "required"}
        sample: {description: "Name of the sample.", category: "required"}
        outputPrefix: {description: "The prefix of the output filename before the .bam extension.", category: "advanced"}
        referenceMMI: {description: "MMI file for the reference.", category: "required"}
        queryFile: {description: "BAM file with reads to align against the reference.", category: "required"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}
        unmapped: {description: "Include unmapped reads in the output.", category: "common"}

        sortThreads: {description: "Extra sorting threads used for samtools sort", category: "advanced"}
        sortMemoryGb: {description: "Amount of memory set for sorting", category: "advanced"}
        cores: {description: "The number of cores to be used.", category: "advanced"}
        memory: {description: "The amount of memory available to the job.", category: "advanced"}

        # output
        outputAlignmentFile: {description: "Mapped bam file."}
        outputIndexFile: {description: "Bam index file."}
    }
}
