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
        mkdir -p "$(dirname ~{outputTable})"
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

    parameter_meta {
        inputBams: {
            description: "The input BAM files.",
            category: "required"
        }
        gtfFile: {
            description: "A GTF/GFF file containing the features of interest.",
            category: "required"
        }
        outputTable: {
            description: "The path to which the output table should be written.",
            category: "common"
        }
        format: {
            description: "Equivalent to the -f option of htseq-count.",
            category: "advanced"
        }
        order: {
            description: "Equivalent to the -r option of htseq-count.",
            category: "advanced"
        }
        stranded: {
            description: "Equivalent to the -s option of htseq-count.",
            category: "common"
        }
        featureType: {
            description: "Equivalent to the --type option of htseq-count.",
            category: "advanced"
        }
        idattr: {
            description: "Equivalent to the --idattr option of htseq-count.",
            category: "advanced"
        }
        additionalAttributes: {
            description: "Equivalent to the --additional-attr option of htseq-count.",
            category: "advanced"
        }
        memory: {
            description: "The amount of memory the job requires in GB.",
            category: "advanced"
        }
        dockerImage: {
            description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
            category: "advanced"
        }
    }
}
