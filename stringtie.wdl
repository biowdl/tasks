task Stringtie {
    String? preCommand
    File alignedReads
    File? referenceGFF
    Int? threads
    String assembledTranscriptsFile

    command {
        set -e -o pipefail
        ${preCommand}
        stringtie \
        ${"-p " + threads} \
        ${"-G " + referenceGFF} \
        ${alignedReads} \
        > ${assembledTranscriptsFile}
    }

    output {
        File assembledTranscripts = assembledTranscriptsFile
    }

    runtime {
        threads: threads
    }
}