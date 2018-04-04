task BwaMem {
    String? preCommand
    File inputR1
    File? inputR2
    String referenceFasta
    String outputPath
    String? readgroup

    Int? threads
    Int? memory

    command {
        set -e -o pipefail
        mkdir -p $(dirname ${outputPath})
        ${preCommand}
        bwa mem ${"-t " + threads} \
        ${"-R '" + readgroup + "'"} \
        ${referenceFasta} ${inputR1} ${inputR2} | samtools sort --output-fmt BAM - > ${outputPath}
    }

    output {
        File bamFile = outputPath
    }
    runtime{
        threads: if defined(threads) then threads else 1
        memory: if defined(memory) then memory else 8
    }
}
