task sample {
    File sequenceFile
    String? outFilePath = "subsampledReads.fq.gz"
    String? preCommand
    Int? seed
    Boolean? twoPassMode
    Float? fraction
    Int? number
    Boolean? zip = true

    command {
    set -e -o pipefail
    ${'mkdir -p $(dirname ' + outFilePath + ')'}
    ${preCommand}
    seqtk sample \
    ${"-s " + seed} \
    ${true="-2 " false="" twoPassMode} \
    ${sequenceFile} \
    ${number} ${fraction} \
    ${true="| gzip" false="" zip} \
    ${"> " + outFilePath}
    }
    output {
        File subsampledReads= select_first([outFilePath])
    }
}