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

task Bam2Fasta {
    input {
        Array[File]+ bam
        Array[File]+ bamIndex
        String outputPrefix
        Int compressionLevel = 1
        Boolean splitByBarcode = false

        String? seqIdPrefix

        String memory = "2G"
        Int timeMinutes = 15
        String dockerImage = "quay.io/biocontainers/bam2fastx:1.3.1--hf05d43a_1"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPrefix})"

        # Localise the bam and pbi files so they are next to each other in the
        # current folder.
        bamFiles=""
        for bamFile in ~{sep=" " bam}
        do
            ln $bamFile .
            bamFiles=$bamFiles" $(basename $bamFile)"
        done

        for index in ~{sep=" " bamIndex}
        do
            ln $index .
        done

        bam2fasta \
        --output ~{outputPrefix} \
        -c ~{compressionLevel} \
        ~{true="--split-barcodes" false="" splitByBarcode} \
        ~{"--seqid-prefix " + seqIdPrefix} \
        $bamFiles
    }

    output {
        File fastaFile = outputPrefix + ".fasta.gz"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        bam: {description: "The input pacbio bam file(s).", category: "required"}
        bamIndex: {description: "The .pbi index for the input file(s).", category: "required"}
        outputPrefix: {description: "Output directory path + output file prefix.", category: "required"}
        compressionLevel: {description: "Gzip compression level [1-9]", category: "advanced"}
        splitByBarcode: {description: "Split output into multiple fasta files, by barcode pairs.", category: "advanced"}
        seqIdPrefix: {description: "Prefix for sequence IDs in headers.", category: "advanced"}
        memory: {description: "The amount of memory available to the job.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        fastaFile: {description: "The fasta output file."}
    }
}

task Bam2Fastq {
    input {
        Array[File]+ bam
        Array[File]+ bamIndex
        String outputPrefix
        Int compressionLevel = 1
        Boolean splitByBarcode = false

        String? seqIdPrefix

        String memory = "2G"
        Int timeMinutes = 15
        String dockerImage = "quay.io/biocontainers/bam2fastx:1.3.1--hf05d43a_1"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPrefix})"

        # Localise the bam and pbi files so they are next to each other in the
        # current folder.
        bamFiles=""
        for bamFile in ~{sep=" " bam}
        do
            ln $bamFile .
            bamFiles=$bamFiles" $(basename $bamFile)"
        done

        for index in ~{sep=" " bamIndex}
        do
            ln $index .
        done

        bam2fastq \
        --output ~{outputPrefix} \
        -c ~{compressionLevel} \
        ~{true="--split-barcodes" false="" splitByBarcode} \
        ~{"--seqid-prefix " + seqIdPrefix} \
        $bamFiles
    }

    output {
        File fastqFile = outputPrefix + ".fastq.gz"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        bam: {description: "The input pacbio bam file(s).", category: "required"}
        bamIndex: {description: "The .pbi index for the input file(s).", category: "required"}
        outputPrefix: {description: "Output directory path + output file prefix.", category: "required"}
        compressionLevel: {description: "Gzip compression level [1-9]", category: "advanced"}
        splitByBarcode: {description: "Split output into multiple fastq files, by barcode pairs.", category: "advanced"}
        seqIdPrefix: {description: "Prefix for sequence IDs in headers.", category: "advanced"}
        memory: {description: "The amount of memory available to the job.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        fastqFile: {description: "The fastq output file."}
    }
}
