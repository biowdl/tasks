version 1.0

task Mem {
    input {
        File read1
        File? read2
        BwaIndex bwaIndex
        String outputPath
        String? readgroup

        Int threads = 2
        String memory = "16G"
        String picardXmx = "4G"
        # A mulled container is needed to have both picard and bwa in one container.
        # This container contains: picard (2.18.7), bwa (0.7.17-r1188)
        String dockerImage = "quay.io/biocontainers/mulled-v2-002f51ea92721407ef440b921fb5940f424be842:43ec6124f9f4f875515f9548733b8b4e5fed9aa6-0"
    }

    command {
        set -e -o pipefail
        mkdir -p $(dirname ~{outputPath})
        bwa mem \
        ~{"-t " + threads} \
        ~{"-R '" + readgroup}~{true="'" false="" defined(readgroup)} \
        ~{bwaIndex.fastaFile} \
        ~{read1} \
        ~{read2} \
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
        docker: dockerImage
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

        Int threads = 2
        Int sortThreads = 2
        String memory = "10G"
        String dockerImage = "biocontainers/bwakit:v0.7.15_cv1"
    }

    command {
        set -e
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
          -@ ~{sortThreads} \
          -m1G \
          - \
          -o ~{outputPrefix}.aln.bam
        samtools index ~{outputPrefix}.aln.bam ~{outputPrefix}.aln.bai
    }

    output {
        File outputBam = outputPrefix + ".aln.bam"
        File outputBamIndex = outputPrefix + ".aln.bai"
    }

    runtime {
        cpu: threads + sortThreads
        memory: memory
        docker: dockerImage
    }

    parameter_meta {
        read1: {
            description: "The first-end fastq file.",
            category: "required"
        }
        read2: {
            description: "The second-end fastq file.",
            category: "common"
        }
        bwaIndex: {
            description: "The BWA index, including a .alt file.",
            category: "required"
        }
        outputPrefix: {
            description: "The prefix of the output files, including any parent directories.",
            category: "required"
        }
        readgroup: {
            description: "A readgroup identifier.",
            category: "common"
        }
        sixtyFour: {
            description: "Whether or not the index uses the '.64' suffixes.",
            category: "common"
        }
        threads: {
            description: "The number of threads to use for alignment.",
            category: "advanced"
        }
        sortThreads: {
            description: "The number of threads to use for sorting.",
            category: "advanced"
        }
        memory: {
            description: "The amount of memory this job will use.",
            category: "advanced"
        }
        dockerImage: {
            description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
            category: "advanced"
        }

        outputBam: "The produced BAM file."
        outputBamIndex: "The index of the produced BAM file."
    }
}

struct BwaIndex {
    File fastaFile
    Array[File] indexFiles
}
