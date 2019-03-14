version 1.0

import "common.wdl"

task Stringtie {
    input {
        String? preCommand
        IndexedBamFile bamFile
        File? referenceGtf
        Int threads = 1
        String assembledTranscriptsFile
        Boolean skipNovelTranscripts = true
        Boolean? firstStranded
        Boolean? secondStranded
        String? geneAbundanceFile
    }

    command {
        set -e -o pipefail
        mkdir -p $(dirname ~{assembledTranscriptsFile})
        ~{preCommand}
        stringtie \
        ~{"-p " + threads} \
        ~{"-G " + referenceGtf} \
        ~{true="-e" false="" skipNovelTranscripts} \
        ~{true="--rf" false="" firstStranded} \
        ~{true="--fr" false="" secondStranded} \
        -o ~{assembledTranscriptsFile} \
        ~{"-A " + geneAbundanceFile} \
        ~{bamFile.file}
    }

    output {
        File assembledTranscripts = assembledTranscriptsFile
        File? geneAbundance = geneAbundanceFile
    }

    runtime {
        cpu: threads
    }
}