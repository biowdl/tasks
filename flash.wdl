version 1.0

task Flash {
    input {
        String? preCommand
        File inputR1
        File inputR2
        String outdirPath
        String outPrefix = "flash"
        Int? minOverlap
        Int? maxOverlap
        Boolean compress = true
        Int threads = 2
        Int memory = 2
    }

    command {
        set -e -o pipefail
        mkdir -p ~{outdirPath}
        ~{preCommand}
        flash \
        ~{"--threads=" + threads} \
        ~{"--output-directory=" + outdirPath} \
        ~{"--output-prefix=" + outPrefix} \
        ~{true="--compress " false="" compress} \
        ~{"--min-overlap=" + minOverlap} \
        ~{"--max-overlap=" + maxOverlap} \
        ~{inputR1} ~{inputR2}
    }

    output {
        File extendedFrags = outdirPath + "/" + outPrefix + ".extendedFrags.fastq.gz"
        File notCombined1 = outdirPath + "/" + outPrefix + ".notCombined_1.fastq.gz"
        File notCombined2 = outdirPath + "/" + outPrefix + ".notCombined_2.fastq.gz"
        File hist = outdirPath + "/" + outPrefix + ".hist"
        File histogram = outdirPath + "/" + outPrefix + ".histogram"
    }

    runtime {
        cpu: threads
        memory: memory
    }

}