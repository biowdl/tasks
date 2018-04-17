task Stringtie {
    String? preCommand
    File alignedReads
    File? referenceGtf
    Int? threads
    String assembledTranscriptsFile
    Boolean? firstStranded
    Boolean? secondStranded
    String? geneAbundanceFile

    command {
        set -e -o pipefail
        ${preCommand}
        stringtie \
        ${"-p " + threads} \
        ${"-G " + referenceGtf} \
        ${true="--rf" false="" firstStranded} \
        ${true="fr" false="" secondStranded} \
        -o ${assembledTranscriptsFile} \
        ${"-A " + geneAbundanceFile} \
        ${alignedReads} \

    }

    output {
        File assembledTranscripts = assembledTranscriptsFile
        File? geneAbundance = geneAbundanceFile
    }

    runtime {
        cpu: select_first([threads, 1])
    }
}