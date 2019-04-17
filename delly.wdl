version 1.0

task CallSV {
    input {
        IndexedBamFile bamFile
        Reference reference
        String outputPath        
    }
    

    command <<<
        delly call \
        -o ~{outputPath} \
        -g ~{reference.fasta} \
        ~{bamFile.file}
    >>>

    output {
        File dellyVcf = "~{outputPath}" 
    }
    
    runtime {
        docker: "quay.io/biocontainers/delly:0.8.1--h4037b6b_1"
    }

}
