version 1.0

task HTSeqCount {
    input {
        Array[File]+ inputBams
        File gtfFile
        String outputTable = "output.tsv"
        String format = "bam"
        String order = "pos"
        String stranded = "no"
        String? featureType
        String? idattr
        Array[String] additionalAttributes = []

        String memory = "40G"
        String dockerImage = "quay.io/biocontainers/htseq:0.9.1--py36h7eb728f_2"
    }

    command {
        set -e
        mkdir -p $(dirname ~{outputTable})
        htseq-count \
        -f ~{format} \
        -r ~{order} \
        -s ~{stranded} \
        ~{"--type " + featureType} \
        ~{"--idattr " + idattr} \
        ~{true="--additional-attr " false="" length(additionalAttributes) > 0 }~{sep=" --additional-attr " additionalAttributes} \
        ~{sep=" " inputBams} \
        ~{gtfFile} \
        > ~{outputTable}
    }

    output {
        File counts = outputTable
    }

    runtime {
        memory: memory
        docker: dockerImage
    }
}
