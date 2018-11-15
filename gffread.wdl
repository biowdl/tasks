version 1.0

task GffRead {
    input {
        String? preCommand
        File inputGff
        File genomicSequence
        File? genomicIndex  # Optional. GFFRead can create this by itself.
        String? exonsFastaPath
        String? CDSFastaPath
        String? proteinFastaPath
        String? filteredGffPath
        Boolean outputGtfFormat = false
    }

    command {
        set -e -o pipefail
        ~{preCommand}
        ~{"mkdir -p $(dirname " + exonsFastaPath +")"}
        ~{"mkdir -p $(dirname " + CDSFastaPath +")"}
        ~{"mkdir -p $(dirname " + proteinFastaPath +")"}
        ~{"mkdir -p $(dirname " + filteredGffPath +")"}
        gffread \
        ~{inputGff} \
        -g ~{genomicSequence} \
        ${"-w " + exonsFastaPath} \
        ${"-x " + CDSFastaPath} \
        ${"-y " + proteinFastaPath} \
        ${"-o " + filteredGffPath} \
        ${true="-T " false="" outputGtfFormat}
    }

    output {
        File? exonsFasta = exonsFastaPath
        File? CDSFasta = CDSFastaPath
        File? proteinFasta = proteinFastaPath
        File? filteredGff = filteredGffPath
    }
}