version 1.0

import "common.wdl"

task Merge {
    input{
        String dockerImage = "quay.io/biocontainers/survivor:1.0.6--h6bb024c_0"
        Array[File] filePaths
        Int breakpointDistance = 1000
        Int suppVecs = 2
        Int svType = 1
        Int strandType = 1
        Int distanceBySvSize = 0
        Int minSize = 30
        String sample
        String outputPath
        Int memory = 128
    }

    command { 
        set -e
        mkdir -p "$(dirname ~{outputPath})"
        echo '~{sep="\n" filePaths}' > fileList
        SURVIVOR merge \
        fileList \
        ~{breakpointDistance} \
        ~{suppVecs} \
        ~{svType} \
        ~{strandType} \
        ~{distanceBySvSize} \
        ~{minSize} \
        ~{outputPath}
    } 

    output {
        File mergedVcf = outputPath
    }
    
    runtime {
        docker: dockerImage
        memory: memory
    }
}
