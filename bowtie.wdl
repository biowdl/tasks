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

task Bowtie {
    input {
        Array[File]+ readsUpstream
        Array[File] readsDownstream = []
        String outputPath = "mapped.bam"
        Array[File]+ indexFiles
        Int? seedmms
        Int? seedlen
        Int? k
        Boolean best = false
        Boolean strata = false
        Boolean allowContain = false
        String? samRG

        Int threads = 1
        Int timeMinutes = 1 + ceil(size(flatten([readsUpstream, readsDownstream]), "G") * 300 / threads)
        String memory = "10G"
        String picardXmx = "4G"
        # Image contains bowtie=1.2.2 and picard=2.9.2
        String dockerImage = "quay.io/biocontainers/mulled-v2-bfe71839265127576d3cd749c056e7b168308d56:1d8bec77b352cdcf3e9ff3d20af238b33ed96eae-0"
    }

    # Assume fastq input with -q flag.
    # The output always needs to be SAM as it is piped into Picard SortSam
    # Hence, the --sam flag is used.

    command {
        set -e -o pipefail
        mkdir -p "$(dirname ~{outputPath})"
        bowtie -q \
        --sam \
        ~{"--seedmms " +  seedmms} \
        ~{"--seedlen " + seedlen} \
        ~{"-k " + k} \
        ~{true="--best" false="" best} \
        ~{true="--strata" false="" strata} \
        ~{true="--allow-contain" false="" allowContain} \
        ~{"--threads " + threads} \
        ~{"--sam-RG '" + samRG}~{true="'" false="" defined(samRG)} \
        ~{sub(indexFiles[0], "(\.rev)?\.[0-9]\.ebwt$", "")} \
        ~{true="-1" false="" length(readsDownstream) > 0} ~{sep="," readsUpstream} \
        ~{true="-2" false="" length(readsDownstream) > 0} ~{sep="," readsDownstream} \
        | picard -Xmx~{picardXmx} SortSam \
        INPUT=/dev/stdin \
        OUTPUT=~{outputPath} \
        SORT_ORDER=coordinate \
        CREATE_INDEX=true
    }

    output {
        File outputBam = outputPath
        File outputBamIndex = sub(outputPath, "\.bam$", ".bai")
    }

    runtime {
        cpu: threads
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        readsUpstream: {description: "The first-/single-end fastq files.", category: "required"}
        readsDownstream: {description: "The second-end fastq files.", category: "common"}
        outputPath: {description: "The location the output BAM file should be written to.", category: "common"}
        indexFiles: {description: "The index files for bowtie.", category: "required"}
        seedmms: {description: "Equivalent to bowtie's `--seedmms` option.", category: "advanced"}
        seedlen: {description: "Equivalent to bowtie's `--seedlen` option.", category: "advanced"}
        k: {description: "Equivalent to bowtie's `-k` option.", category: "advanced"}
        best: {description: "Equivalent to bowtie's `--best` flag.", category: "advanced"}
        strata: {description: "Equivalent to bowtie's `--strata` flag.", category: "advanced"}
        allowContain: {description: "Equivalent to bowtie's `--allow-contain` flag.", category: "advanced"}
        samRG: {description: "Equivalent to bowtie's `--sam-RG` option.", category: "advanced"}

        picardXmx: {description: "The maximum memory available to the picard (used for sorting the output). Should be lower than `memory` to accommodate JVM overhead and bowtie's memory usage.",
                  category: "advanced"}
        threads: {description: "The number of threads to use.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}
