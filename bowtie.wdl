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
        Int memory = 8
        Int picardMemory = 4
        # Image contains bowtie=1.2.2 and picard=2.9.2
        String dockerImage = "quay.io/biocontainers/mulled-v2-bfe71839265127576d3cd749c056e7b168308d56:1d8bec77b352cdcf3e9ff3d20af238b33ed96eae-0"
    }

    # Assume fastq input with -q flag.
    # The output always needs to be SAM as it is piped into Picard SortSam
    # Hence, the --sam flag is used.

    command {
        set -e -o pipefail
        mkdir -p $(dirname ~{outputPath})
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
        | picard -Xmx~{picardMemory}G SortSam \
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
        memory: memory + picardMemory + picardMemory
        docker: dockerImage
    }
}

struct BowtieIndex {
    File fasta
    Array[File] indexFiles
}