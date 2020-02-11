version 1.0

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

        String memory = "24G"
        String dockerImage = "quay.io/biocontainers/survivor:1.0.6--h6bb024c_0"
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
        memory: memory
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        filePaths: {description: "An array of VCF files (predictions) to be merged by SURVIVOR", category: "advanced"}
        breakpointDistance: {description: "The distance between pairwise breakpoints between SVs", category: "advanced"}
        suppVecs: {description: "The minimum number of SV callers to support the merging", category: "advanced"}
        svType: {description: "A boolean to include the type SV to be merged", category: "advanced"}
        strandType: {description: "A boolean to include strand type of an SV to be merged", category: "advanced"}
        distanceBySvSize: {description: "A boolean to predict the pairwise distance between the SVs based on their size", category: "advanced"}
        minSize: {description: "The mimimum size of SV to be merged", category: "advanced"}
        sample: {description: "The name of the sample", category: "required"}
        outputPath: {description: "The location the output VCF file should be written.", category: "common"}
        memory: {description: "The memory required to run the programs", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}
    }
}
