task HTSeqCount {
    String? preCommand
    Array[File] alignmentFiles
    File gffFile
    String outputTable
    String? format
    String? order
    String? stranded

    command {
        set -e -o pipefail
        ${preCommand}
        htseq-count \
        -f ${default="bam" format} \
        -r ${default="pos" order} \
        -s ${default="no" stranded} \
        ${sep=" " alignmentFiles} \
        ${gffFile} \
        > ${outputTable}
    }

    output {
        File counts = outputTable
    }

    runtime {
        memory: 3
    }
}