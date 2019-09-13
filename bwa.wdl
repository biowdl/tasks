version 1.0

task Mem {
    input {
        File read1
        File? read2
        BwaIndex bwaIndex
        String outputPath
        String? readgroup

        Int threads = 2
        Int memory = 8
        Int picardMemory = 4
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

struct BwaIndex {
    File fastaFile
    Array[File] indexFiles
}
