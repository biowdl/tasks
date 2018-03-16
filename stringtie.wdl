task Stringtie {
    String? preCommand
    File alignedReads
    File? referenceGFF
    Int? threads
    String assembledTranscriptsFile
    Boolean? firstStranded
    Boolean? secondStranded

    command {
        set -e -o pipefail
        ${preCommand}
        stringtie \
        ${"-p " + threads} \
        ${"-G " + referenceGFF} \
        ${true="--rf" false="" firstStranded} \
        ${true="fr" false="" secondStranded} \
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