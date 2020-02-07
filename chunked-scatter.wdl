version 1.0

task ChunkedScatter {
    input {
        File inputFile
        String prefix = "./scatter"
        Int? chunkSize
        Int? overlap
        Int? minimumBasesPerFile

        String dockerImage = "quay.io/biocontainers/chunked-scatter:0.1.0--py_0"
    }

    command {
        set -e
        mkdir -p ~{prefix}
        chunked-scatter \
        -p ~{prefix} \
        -i ~{inputFile} \
        ~{"-c " + chunkSize} \
        ~{"-o " + overlap} \
        ~{"-m " + minimumBasesPerFile}
    }

    output {
        Array[File] scatters = glob(prefix + "*.bed")
    }

    runtime {
        memory: "4G"
        docker: dockerImage
    }

    parameter_meta {
        inputFile: {description: "Either a bed file describing regiosn of intrest or a sequence dictionary.", category: "required"}
        prefix: {description: "The prefix for the output files.", category: "advanced"}
        chunkSize: {description: "Equivalent to chunked-scatter's `-c` option.", category: "advanced"}
        overlap: {description: "Equivalent to chunked-scatter's `-o` option.", category: "advanced"}
        minimumBasesPerFile: {description: "Equivalent to chunked-scatter's `-m` option.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}