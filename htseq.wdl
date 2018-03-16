task HTSeqCount {
    Array[File] alignmentFiles
    File gffFile

    String? format
    String? order
    String? stranded

    String? preCommand

    command {
        set -e -o pipefail
        ${preCommand}
        htseq-count \
        -f ${default="bam" format} \
        -r ${default="pos" order} \
        -s ${default="no" stranded} \
        ${sep=" " alignmentFiles} \
        ${gffFile}
    }

    output {
        File counts = stdout()
    }
}