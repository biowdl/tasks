version 1.0

task Awk {
    input {
        File in
        String outputPrefix
        String awk
        String? sep

        Int threads = 1
        Int timeMinutes = 10 + size(in, "GiB")
        # Contains bwa 0.7.17 bwakit 0.7.17.dev1 and samtools 1.10.
        String dockerImage = "quay.io/biocontainers/samtools:1.21--h96c455f_1"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPrefix})"

        cat ~{in} | \
            awk \
                ~{"-F " + sep} \
                ~{awk} \
                > ~{outputPrefix}
    }

    output {
        File out = outputPrefix
    }

    runtime {
        cpu: threads
        memory: "1GiB"
        time_minutes: timeMinutes
        docker: dockerImage
    }
}
