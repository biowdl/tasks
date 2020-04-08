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
        String memory = "20G"
        String picardXmx = "4G"
        Int timeMinutes = 1 + ceil(size([read1, read2], "G") * 200 / threads)
        # A mulled container is needed to have both picard and bwa in one container.
        # This container contains: picard (2.18.7), bwa (0.7.17-r1188)
        String dockerImage = "quay.io/biocontainers/mulled-v2-002f51ea92721407ef440b921fb5940f424be842:43ec6124f9f4f875515f9548733b8b4e5fed9aa6-0"
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
        | picard -Xmx~{picardXmx} -XX:ParallelGCThreads=1 SortSam \
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
        read1: {description: "The first or single end fastq file.", category: "required"}
        read2: {description: "The second end fastq file.", category: "common"}
        bwaIndex: {description: "The BWA index files.", category: "required"}
        outputPath: {description: "The location the output BAM file should be written to.", category: "required"}
        readgroup: {description: "The readgroup to be assigned to the reads. See BWA mem's `-R` option.", category: "common"}

        threads: {description: "The number of threads to use.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        picardXmx: {description: "The maximum memory available to picard SortSam. Should be lower than `memory` to accommodate JVM overhead and BWA mem's memory usage.",
                  category: "advanced"}
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
        # Samtools uses *additional* threads. So by default this option should
        # not be used.
        Int? sortThreads
        # Compression uses zlib. Higher than level 2 causes enormous slowdowns.
        # GATK/Picard default is level 2.
        String sortMemoryPerThread = "4G"
        Int compressionLevel = 1
        String memory = "20G"
        Int timeMinutes = 1 + ceil(size([read1, read2], "G") * 220 / threads)
        String dockerImage = "biocontainers/bwakit:v0.7.15_cv1"
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
        k8 /opt/conda/bin/bwa-postalt.js \
          -p ~{outputPrefix}.hla \
          ~{bwaIndex.fastaFile}~{true=".64.alt" false=".alt" sixtyFour} | \
        samtools sort \
          ~{"-@ " + sortThreads} \
          -m ~{sortMemoryPerThread} \
          -l ~{compressionLevel} \
          - \
          -o ~{outputPrefix}.aln.bam
        samtools index ~{outputPrefix}.aln.bam ~{outputPrefix}.aln.bai
    }

    output {
        File outputBam = outputPrefix + ".aln.bam"
        File outputBamIndex = outputPrefix + ".aln.bai"
    }

    runtime {
        cpu: threads + 1  # One thread for bwa-postalt + samtools.
        memory: memory
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
        sortThreads: {description: "The number of additional threads to use for sorting.", category: "advanced"}
        sortMemoryPerThread: {description: "The amount of memory for each sorting thread.", category: "advanced"}
        compressionLevel: {description: "The compression level of the output BAM.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}

        # outputs
        outputBam: "The produced BAM file."
        outputBamIndex: "The index of the produced BAM file."
    }
}

struct BwaIndex {
    File fastaFile
    Array[File] indexFiles
}
