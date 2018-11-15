version 1.0

task gffread {
    input {
        String? preCommand
        File inputGff
        File genomicSequences
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
        -g ~{genomicSequences} \
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