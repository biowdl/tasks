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

task Mem {
    input {
        File read1
        File? read2
        BwaIndex bwaIndex
        String outputPath
        String? readgroup

        Int threads = 4
        Int sortThreads = 1
        Int sortMemoryPerThreadGb = 2
        Int compressionLevel = 1
        # BWA needs slightly more memory than the size of the index files (~10%). Add a margin for safety here.
        Int memoryGb = 1 + ceil(size(bwaIndex.indexFiles, "G") * 1.2) + sortMemoryPerThreadGb * sortThreads
        Int timeMinutes = 1 + ceil(size([read1, read2], "G") * 200 / threads)
        # This container contains: samtools (1.10), bwa (0.7.17-r1188)
        String dockerImage = "quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:eabfac3657eda5818bae4090db989e3d41b01542-0"
    }

    command {
        set -e -o pipefail
        mkdir -p "$(dirname ~{outputPath})"
        bwa mem \
        ~{"-t " + threads} \
        ~{"-R '" + readgroup}~{true="'" false="" defined(readgroup)} \
        ~{bwaIndex.fastaFile} \
        ~{read1} \
        ~{read2} \
        | samtools sort \
        ~{"-@ " + sortThreads} \
        -m ~{sortMemoryPerThreadGb}G \
        -l ~{compressionLevel} \
        - \
        -o ~{outputPath}
    }

    output {
        File outputBam = outputPath
    }

    runtime {
        cpu: threads
        memory: "~{memoryGb}G"
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        read1: {description: "The first or single end fastq file.", category: "required"}
        read2: {description: "The second end fastq file.", category: "common"}
        bwaIndex: {description: "The BWA index files.", category: "required"}
        outputPath: {description: "The location the output BAM file should be written to.", category: "required"}
        readgroup: {description: "The readgroup to be assigned to the reads. See BWA mem's `-R` option.", category: "common"}

        threads: {description: "The number of threads to use.", category: "advanced"}
        memoryGb: {description: "The amount of memory this job will use in gigabytes.", category: "advanced"}
        sortThreads: {description: "The number of threads to use for sorting.", category: "advanced"}
        sortMemoryPerThreadGb: {description: "The amount of memory for each sorting thread in gigabytes.", category: "advanced"}
        compressionLevel: {description: "The compression level of the output BAM.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}

task Kit {
    input {
        File read1
        File? read2
        BwaIndex bwaIndex
        String outputPrefix
        String? readgroup
        Boolean sixtyFour = false

        Int threads = 4
        Int sortThreads = 1
        Int sortMemoryPerThreadGb = 2
        Int compressionLevel = 1
        # BWA needs slightly more memory than the size of the index files (~10%). Add a margin for safety here.
        Int memoryGb = 1 + ceil(size(bwaIndex.indexFiles, "G") * 1.2) + sortMemoryPerThreadGb * sortThreads
        Int timeMinutes = 1 + ceil(size([read1, read2], "G") * 220 / threads)
        # Contains bwa 0.7.17 bwakit 0.7.17.dev1 and samtools 1.10
        String dockerImage = "quay.io/biocontainers/mulled-v2-ad317f19f5881324e963f6a6d464d696a2825ab6:c59b7a73c87a9fe81737d5d628e10a3b5807f453-0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPrefix})"
        bwa mem \
          -t ~{threads} \
          ~{"-R '" + readgroup}~{true="'" false="" defined(readgroup)} \
          ~{bwaIndex.fastaFile} \
          ~{read1} \
          ~{read2} \
          2> ~{outputPrefix}.log.bwamem | \
        bwa-postalt.js \
          -p ~{outputPrefix}.hla \
          ~{bwaIndex.fastaFile}~{true=".64.alt" false=".alt" sixtyFour} | \
        samtools sort \
          ~{"-@ " + sortThreads} \
          -m ~{sortMemoryPerThreadGb}G \
          -l ~{compressionLevel} \
          - \
          -o ~{outputPrefix}.aln.bam
    }

    output {
        File outputBam = outputPrefix + ".aln.bam"
    }

    runtime {
        # One extra thread for bwa-postalt + samtools is not needed.
        # These only use 5-10% of compute power and not always simultaneously.
        cpu: threads  
        memory: "~{memoryGb}G"
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        read1: {description: "The first-end fastq file.", category: "required"}
        read2: {description: "The second-end fastq file.", category: "common"}
        bwaIndex: {description: "The BWA index, including a .alt file.", category: "required"}
        outputPrefix: {description: "The prefix of the output files, including any parent directories.", category: "required"}
        readgroup: {description: "A readgroup identifier.", category: "common"}
        sixtyFour: {description: "Whether or not the index uses the '.64' suffixes.", category: "common"}
        threads: {description: "The number of threads to use for alignment.", category: "advanced"}
        memoryGb: {description: "The amount of memory this job will use in gigabytes.", category: "advanced"}
        sortThreads: {description: "The number of threads to use for sorting.", category: "advanced"}
        sortMemoryPerThreadGb: {description: "The amount of memory for each sorting thread in gigabytes.", category: "advanced"}
        compressionLevel: {description: "The compression level of the output BAM.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}

        # outputs
        outputBam: "The produced BAM file."
    }
}

struct BwaIndex {
    File fastaFile
    Array[File] indexFiles
}
