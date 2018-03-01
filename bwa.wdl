task BwaMem {
    File inputR1
    File? inputR2
    String referenceFasta
    String outputPath
    String? readgroup

    command {
        set -e -o pipefail
        mkdir -p $(dirname ${outputPath})
        bwa mem ${"-R '" + readgroup + "'"} \
        ${referenceFasta} ${inputR1} ${inputR2} | samtools sort --output-fmt BAM - > ${outputPath}
    }

    output {
        File bamFile = outputPath
    }
}
