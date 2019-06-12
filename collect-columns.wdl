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

        String dockerTag = "0.2.0--py_1"
    }

    command {
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

    runtime {
        memory: 4 + ceil(0.5* length(inputTables))
        docker: "quay.io/biocontainers/collect-columns:" + dockerTag
    }
}