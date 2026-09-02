version 1.0

task Awk {
    input {
        File in
        String outputPrefix
        String awk
        String? ofs
        String? sep

        Int threads = 1
        Int timeMinutes = 10 + ceil(size(in, "GiB"))
        # Contains bwa 0.7.17 bwakit 0.7.17.dev1 and samtools 1.10.
        String dockerImage = "quay.io/biocontainers/samtools:1.21--h96c455f_1"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPrefix})"

        cat ~{in} | \
            awk ~{"-F " + sep} ~{"-v OFS=" + ofs} \
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

task Sed {
    input {
        File in
        String outputPrefix
        Boolean extendedRegex = true
        # You are responsible for quoting it yourself!!!!!
        String sed

        Int threads = 1
        Int timeMinutes = 10 + ceil(size(in, "GiB"))
        # Contains bwa 0.7.17 bwakit 0.7.17.dev1 and samtools 1.10.
        String dockerImage = "quay.io/biocontainers/samtools:1.21--h96c455f_1"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPrefix})"

        cat ~{in} | \
            sed ~{true="-E" false="" extendedRegex} \
                ~{sed} \
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
