version 1.0

import "common.wdl"

task CallSV {
    input {
        IndexedBamFile bamFile
        Reference reference
        String outputPath        
        Int mem = 15
    }
    

    command <<<
        set -e
        mkdir -p $(dirname ~{outputPath})
        delly call \
        -o ~{outputPath} \
        -g ~{reference.fasta} \
        ~{bamFile.file}
    >>>

    output {
        File dellyBcf = "~{outputPath}" 
    }
    
    runtime {
        docker: "quay.io/biocontainers/delly:0.8.1--h4037b6b_1"
        memory: mem
    }

}
