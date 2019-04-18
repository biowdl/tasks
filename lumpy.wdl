version 1.0

import "common.wdl"

task CallSV {
    input {
        IndexedBamFile bamFile
        Reference reference
        String outputPath    
    }   
    

    command <<< 
        lumpyexpress \
        -B ~{bamFile.file} \
        -o ~{outputPath}
    >>> 

    output {
        File lumpyVcf = "~{outputPath}"
    }   
        
    runtime {
        docker: "quay.io/biocontainers/lumpy-sv:0.3.0--h0b85cd1_2"
    }   

}
