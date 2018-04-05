task HTSeqCount {
    String? preCommand
    Array[File] alignmentFiles
    File gtfFile
    String outputTable
    String? format
    String? order
    String? stranded

    Int? memory

    command {
        set -e -o pipefail
        ${preCommand}
        htseq-count \
        -f ${default="bam" format} \
        -r ${default="pos" order} \
        -s ${default="no" stranded} \
        ${sep=" " alignmentFiles} \
        ${gtfFile} \
        > ${outputTable}
    }

    output {
        File counts = outputTable
    }

    runtime {
        memory: select_first([memory, 3])
    }
}