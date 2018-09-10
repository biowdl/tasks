version 1.0

import "common.wdl"

task HTSeqCount {
    input {
        String? preCommand
        Array[IndexedBamFile] inputBamFiles
        File gtfFile
        String outputTable
        String format = "bam"
        String order = "pos"
        String stranded = "no"

        Int memory = 3
    }

    command {
        set -e -o pipefail
        mkdir -p ~{sub(outputTable, basename(outputTable), "")}
        ~{preCommand}
        htseq-count \
        -f ~{format} \
        -r ~{order} \
        -s ~{stranded} \
        ~{sep=" " inputBamFiles.file} \
        ~{gtfFile} \
        > ~{outputTable}
    }

    output {
        File counts = outputTable
    }

    runtime {
        memory: memory
    }
}