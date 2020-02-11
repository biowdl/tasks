version 1.0

import "common.wdl"

task CallSV {
    input {
        File bamFile
        File bamIndex
        File referenceFasta
        File referenceFastaFai
        String outputPath

        String memory = "15G"
        String dockerImage = "quay.io/biocontainers/delly:0.8.1--h4037b6b_1"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPath})"
        delly call \
        -o ~{outputPath} \
        -g ~{referenceFasta} \
        ~{bamFile}
    }

    output {
        File dellyBcf = outputPath
    }

    runtime {
        memory: memory
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        bamFile: {description: "The bam file to process.", category: "required"}
        bamIndex: {description: "The index bam file.", category: "required"}
        referenceFasta: {description: "The reference fasta file also used for mapping.", category: "advanced"}
        referenceFastaFai: {description: "Fasta index (.fai) file of the reference", category: "required" }
        outputPath: {description: "The location the output VCF file should be written.", category: "common"}
        memory: {description: "The memory required to run the programs", category: "common"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}
    }
}
