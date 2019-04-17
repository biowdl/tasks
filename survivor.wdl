version 1.0

task Merge {
    input{
        File filePaths
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
        SURVIVOR merge \
        ~{filePaths} \
        ~{breakpointDistance} \
        ~{suppVecs} \
        ~{svType} \
        ~{strandType} \
        ~{distanceBySvSize} \
        ~{minSize} \
        ~{outputPath}
    >>> 

    output {
     File mergedVcf = "~{outputPath}"
    }
    
    runtime {
        docker: "quay.io/biocontainers/survivor:1.0.6--h6bb024c_0"
    }
}
