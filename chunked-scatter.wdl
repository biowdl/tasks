version 1.0

task ChunkedScatter {
    input {
        File inputFile
        String prefix = "./scatter"
        Int? chunkSize
        Int? overlap
        Int? minimumBasesPerFile

        String dockerImage = "alpine:latest" #TODO
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
        Array[File] scatters = glob(prefix + "_*.bed")
    }

    runtime {
        memory: 4
        docker: dockerImage
    }
}