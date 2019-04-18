version 1.0

import "common.wdl"

task Merge {
    input{
        Array[File] filePaths
        Int breakpointDistance = 1000
        Int suppVecs = 2
        Int svType = 1
        Int strandType = 1
        Int distanceBySvSize = 0
        Int minSize = 30
        String sample
        String outputPath
    }

    command <<< 
        set -e
        echo '~{sep="\n" filePaths}' > fileList
        SURVIVOR merge \
        ~{fileList} \
        ~{breakpointDistance} \
        ~{suppVecs} \
        ~{svType} \
        ~{strandType} \
        ~{distanceBySvSize} \
        ~{minSize} \
        ~{outputPath}
    >>> 

    output {
        File fileList = "~{fileList}"
        File mergedVcf = "~{outputPath}"
    }
    
    runtime {
        docker: "quay.io/biocontainers/survivor:1.0.6--h6bb024c_0"
    }
}
