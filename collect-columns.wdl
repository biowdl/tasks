version 1.0

task CollectColumns {
    input {
        Array[File]+ inputTables
        String outputPath
        Int? featureColumn
        Int? valueColumn
        Int? separator
        Array[String]? sampleNames
        Boolean header = false
        Array[String]? additionalAttributes
        File? referenceGtf
        String? featureAttribute

        String dockerImage = "quay.io/biocontainers/collect-columns:0.2.0--py_1"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPath})"
        collect-columns \
        ~{outputPath} \
        ~{sep=" " inputTables} \
        ~{"-f "  + featureColumn} \
        ~{"-c " + valueColumn} \
        ~{"-s " + separator} \
        ~{true="-n" false="" defined(sampleNames)} ~{sep=" " sampleNames} \
        ~{true="-H" false="" header} \
        ~{true="-a" false="" defined(additionalAttributes)} ~{sep=" " additionalAttributes} \
        ~{"-g " + referenceGtf} \
        ~{"-F " + featureAttribute}
    }

    output {
        File outputTable = outputPath
    }

    Int memoryGb = 4 + ceil(0.5 * length(inputTables))

    runtime {
        memory: "~{memoryGb}G"
        docker: dockerImage
    }

    parameter_meta {
        inputTables: {description: "The tables from which columns should be taken.",
                      category: "required"}
        outputPath: {description: "The path to which the output should be written.",
                     category: "required"}
        featureColumn: {description: "Equivalent to the -f option of collect-columns.",
                        category: "advanced"}
        valueColumn: {description: "Equivalent to the -c option of collect-columns.",
                      category: "advanced"}
        separator: {description: "Equivalent to the -s option of collect-columns.",
                    category: "advanced"}
        sampleNames: {description: "Equivalent to the -n option of collect-columns.",
                      category: "advanced"}
        header: {description: "Equivalent to the -H flag of collect-columns.",
                 category: "advanced"}
        additionalAttributes: {description: "Equivalent to the -a option of collect-columns.",
                               category: "advanced"}
        referenceGtf: {description: "Equivalent to the -g option of collect-columns.",
                       category: "advanced"}
        featureAttribute: {description: "Equivalent to the -F option of collect-columns.",
                           category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}