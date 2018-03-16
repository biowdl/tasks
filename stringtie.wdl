task Stringtie {
    String? preCommand
    File alignedReads
    File? referenceGFF
    Int? threads
    String assembledTranscripts

    command {
        set -e -o pipefail
        ${preCommand}
        stringtie \
        ${"-p " + threads} \
        ${"-G " + referenceGFF} \
        ${alignedReads} \
        > ${assembledTranscripts}
    }

    output {
        File assembledTranscripts = assembledTranscripts
    }

    runtime {
        threads: threads
    }
}