version 1.0

# Copyright (c) 2017 Sequencing Analysis Support Core - Leiden University Medical Center
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

task Extract {
    input {
        File read1
        File? read2
        String bcPattern
        String? bcPattern2
        Boolean threePrime = false
        String read1Output = "umi_extracted_R1.fastq.gz"
        String? read2Output = "umi_extracted_R2.fastq.gz"

        String dockerImage = "quay.io/biocontainers/mulled-v2-509311a44630c01d9cb7d2ac5727725f51ea43af:6089936aca6219b5bb5f54210ac5eb456c7503f2-0"
    }

    command {
        umi_tools extract \
        --stdin ~{read1} \
        ~{"--read2-in " + read2} \
        --bc-pattern ~{bcPattern} \
        ~{"bc-pattern2 " + bcPattern2} \
        ~{true="--3prime" false="" threePrime} \
        --stdout ~{read1Output} \
        ~{if defined(read2) then "--read2-out " + read2Output else ""}
    }

    output {
        File extractedRead1 = read1Output
        File? extractedRead2 = read2Output
    }

    runtime {
        docker: dockerImage
    }

    parameter_meta {
        read1: {description: "The first/single-end fastq file.", category: "required"}
        read2: {description: "The second-end fastq file.", category: "common"}
        bcPattern: {description: "The pattern to be used for UMI extraction. See the umi_tools docs for more information.", category: "required"}
        bcPattern2: {description: "The pattern to be used for UMI extraction in the second-end reads. See the umi_tools docs for more information.", category: "advanced"}
        threePrime: {description: "Whether or not the UMI's are at the reads' 3' end. If false the UMIs are extracted from the 5' end.", category: "advanced"}
        read1Output: {description: "The location to write the first/single-end output fastq file to.", category: "advanced"}
        read2Output: {description: "The location to write the second-end output fastq file to.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}

task Dedup {
    input {
        File inputBam
        File inputBamIndex
        String outputBamPath
        String? statsPrefix
        Boolean paired = true

        String memory = "5G"

        # Use a multi-package-container which includes umi_tools (0.5.5) and samtools (1.9)
        String dockerImage = "quay.io/biocontainers/mulled-v2-509311a44630c01d9cb7d2ac5727725f51ea43af:6089936aca6219b5bb5f54210ac5eb456c7503f2-0"
    }

    String outputBamIndex = sub(outputBamPath, "\.bam$", ".bai")

    command {
        set -e
        mkdir -p "$(dirname ~{outputBamPath})"
        umi_tools dedup \
        --stdin ~{inputBam} \
        --stdout ~{outputBamPath} \
        ~{"--output-stats " + statsPrefix} \
        ~{true="--paired" false="" paired}
        samtools index ~{outputBamPath} ~{outputBamIndex}
    }

    output {
        File deduppedBam = outputBamPath
        File deduppedBamIndex = outputBamIndex
        File? editDistance = select_first([statsPrefix, "stats"]) + "_edit_distance.tsv"
        File? umiStats = select_first([statsPrefix, "stats"]) + "_per_umi.tsv"
        File? positionStats = select_first([statsPrefix, "stats"]) + "_per_umi_per_position.tsv"
    }

    runtime {
        docker: dockerImage
        memory: memory
    }

    parameter_meta {
        inputBam: {description: "The input BAM file.", categrory: "required"}
        inputBamIndex: {description: "The index for the ipnut BAM file.", cateogry: "required"}
        outputBamPath: {description: "The location to write the output BAM file to.", category: "required"}
        statsPrefix: {description: "The prefix for the stats files.", category: "advanced"}
        paired: {description: "Whether or not the data is paired.", category: "common"}
        memory: {description: "The amount of memory required for the task.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}