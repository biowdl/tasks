version 1.0

task GffRead {
    input {
        File inputGff
        File genomicSequence
        File? genomicIndex  # Optional. GFFRead can create this by itself.
        String? exonsFastaPath
        String? CDSFastaPath
        String? proteinFastaPath
        String? filteredGffPath
        Boolean outputGtfFormat = false
        String dockerTag = "0.9.12--0"
    }

    # The mkdirs below are hackish. It should be
    # ~{"mkir -p $(dirname " + somePath + ")"}
    # but this goes wrong. Cromwell will always use ')' even if somepath is not defined.
    # Which leads to crashing.
    command {
        set -e
        ~{"mkdir -p $(dirname " + CDSFastaPath}~{true=")" false="" defined(CDSFastaPath)}
        ~{"mkdir -p $(dirname " + exonsFastaPath}~{true=")" false="" defined(exonsFastaPath)}
        ~{"mkdir -p $(dirname " + proteinFastaPath}~{true=")" false="" defined(proteinFastaPath)}
        ~{"mkdir -p $(dirname " + filteredGffPath}~{true=")" false="" defined(filteredGffPath)}
        gffread \
        ~{inputGff} \
        -g ~{genomicSequence} \
        ~{"-w " + exonsFastaPath} \
        ~{"-x " + CDSFastaPath} \
        ~{"-y " + proteinFastaPath} \
        ~{"-o " + filteredGffPath} \
        ~{true="-T " false="" outputGtfFormat}
    }

    output {
        File? exonsFasta = exonsFastaPath
        File? CDSFasta = CDSFastaPath
        File? proteinFasta = proteinFastaPath
        File? filteredGff = filteredGffPath
    }

    runtime {
        docker: "quay.io/biocontainers/gffread:" + dockerTag
    }
}