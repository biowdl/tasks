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

task Extract {
    input {
        File read1
        File? read2
        String bcPattern
        String? bcPattern2
        String read1Output = "umi_extracted_R1.fastq.gz"
        String? read2Output = "umi_extracted_R2.fastq.gz"
        Boolean threePrime = false

        String memory = "20GiB"
        Int timeMinutes = 1 + ceil(size([read1, read2], "G") * 2)
        String dockerImage = "quay.io/biocontainers/mulled-v2-509311a44630c01d9cb7d2ac5727725f51ea43af:3067b520386698317fd507c413baf7f901666fd4-0"
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
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        read1: {description: "The first/single-end fastq file.", category: "required"}
        read2: {description: "The second-end fastq file.", category: "common"}
        bcPattern: {description: "The pattern to be used for UMI extraction. See the umi_tools docs for more information.", category: "required"}
        bcPattern2: {description: "The pattern to be used for UMI extraction in the second-end reads. See the umi_tools docs for more information.", category: "advanced"}
        read1Output: {description: "The location to write the first/single-end output fastq file to.", category: "advanced"}
        read2Output: {description: "The location to write the second-end output fastq file to.", category: "advanced"}
        threePrime: {description: "Whether or not the UMI's are at the reads' 3' end. If false the UMIs are extracted from the 5' end.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        extractedRead1: {description: "First read with UMI extracted to read name."}
        extractedRead2: {description: "Second read with UMI extracted to read name."}
    }
}

task Dedup {
    input {
        File inputBam
        File inputBamIndex
        String outputBamPath
        String tmpDir = "./umiToolsDedupTmpDir"

        Boolean paired = true

        String? umiSeparator
        String? statsPrefix

        String memory = "25GiB"
        Int timeMinutes = 30 + ceil(size(inputBam, "GiB") * 30)
        String dockerImage = "quay.io/biocontainers/mulled-v2-509311a44630c01d9cb7d2ac5727725f51ea43af:3067b520386698317fd507c413baf7f901666fd4-0"
    }

    String outputBamIndex = sub(outputBamPath, "\.bam$", ".bai")

    command {
        set -e
        mkdir -p "$(dirname ~{outputBamPath})" "~{tmpDir}"
        umi_tools dedup \
        --stdin=~{inputBam} \
        --stdout=~{outputBamPath} \
        ~{"--output-stats " + statsPrefix} \
        ~{"--umi-separator=" + umiSeparator} \
        ~{true="--paired" false="" paired} \
        --temp-dir=~{tmpDir}
        samtools index ~{outputBamPath} ~{outputBamIndex}
    }

    output {
        File deduppedBam = outputBamPath
        File deduppedBamIndex = outputBamIndex
        File? editDistance = "~{statsPrefix}_edit_distance.tsv"
        File? umiStats = "~{statsPrefix}_per_umi.tsv"
        File? positionStats =  "~{statsPrefix}_per_umi_per_position.tsv"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        inputBam: {description: "The input BAM file.", categrory: "required"}
        inputBamIndex: {description: "The index for the ipnut BAM file.", cateogry: "required"}
        outputBamPath: {description: "The location to write the output BAM file to.", category: "required"}
        tmpDir: {description: "Temporary directory.", category: "advanced"}
        paired: {description: "Whether or not the data is paired.", category: "common"}
        umiSeparator: {description: "Seperator used for UMIs in the read names.", category: "advanced"}
        statsPrefix: {description: "The prefix for the stats files.", category: "advanced"}
        memory: {description: "The amount of memory required for the task.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        deduppedBam: {description: "Deduplicated BAM file."}
        deduppedBamIndex: {description: "Index of the deduplicated BAM file."}
        editDistance: {description: "Report of the (binned) average edit distance between the UMIs at each position."}
        umiStats: {description: "UMI-level summary statistics."}
        positionStats: {description: "The counts for unique combinations of UMI and position."}
    }
}
