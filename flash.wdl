version 1.0

import "common.wdl" as common

task Flash {
    input {
        String? preCommand
        FastqPair inputFastq
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
        ~{inputFastq.R1} ~{inputFastq.R2}
    }

    output {
        File extendedFrags = outdirPath + "/" + outPrefix + ".extendedFrags.fastq.gz"
        File notCombined1 = outdirPath + "/" + outPrefix + ".notCombined_1.fastq.gz"
        File notCombined2 = outdirPath + "/" + outPrefix + ".notCombined_2.fastq.gz"
        FastqPair notCombined = object {
          R1: notCombined1,
          R2: notCombined2
        }
        File hist = outdirPath + "/" + outPrefix + ".hist"
        File histogram = outdirPath + "/" + outPrefix + ".histogram"
    }

    runtime {
        cpu: threads
        memory: memory
    }

}