version 1.0

task HTSeqCount {
    input {
        Array[File]+ inputBams
        Array[File]+ inputBamsIndex
        File gtfFile
        String outputTable = "output.tsv"
        String format = "bam"
        String order = "pos"
        String stranded = "no"

        Int memory = 20
        String dockerTag = "0.9.1--py36h7eb728f_2"
    }

    command {
        set -e
        mkdir -p $(dirname ~{outputTable})
        htseq-count \
        -f ~{format} \
        -r ~{order} \
        -s ~{stranded} \
        ~{sep=" " inputBams} \
        ~{gtfFile} \
        > ~{outputTable}
    }

    output {
        File counts = outputTable
    }

    runtime {
        memory: memory
        docker: "quay.io/biocontainers/htseq:" + dockerTag
    }
}