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
}