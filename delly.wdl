version 1.0

import "common.wdl"

task CallSV {
    input {
        File bamFile
        File bamIndex
        File referenceFasta
        File referenceFastaFai
        String outputPath        
        String memory = 15
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
        docker: dockerImage
        memory: memory
    }
}
