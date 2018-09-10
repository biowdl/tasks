version 1.0

import "common.wdl"

task Stringtie {
    input {
        String? preCommand
        IndexedBamFile bam
        File? referenceGtf
        Int threads = 1
        String assembledTranscriptsFile
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
        ~{true="--rf" false="" firstStranded} \
        ~{true="fr" false="" secondStranded} \
        -o ~{assembledTranscriptsFile} \
        ~{"-A " + geneAbundanceFile} \
        ~{bam.file}
    }

    output {
        File assembledTranscripts = assembledTranscriptsFile
        File? geneAbundance = geneAbundanceFile
    }

    runtime {
        cpu: threads
    }
}