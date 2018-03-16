task Stringtie {
    String? preCommand
    File alignedReads
    File? referenceGFF
    Int? threads

    command {
        set -e -o pipefail
        ${preCommand}
        stringtie \
        ${"-p " + threads} \
        ${"-G " + referenceGFF} \
        ${alignedReads}
    }

    output {
        File assembledTranscripts = stdout()
    }

    runtime {
        threads: threads
    }
}