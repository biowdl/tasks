task flash {
    String? preCommand
    File inputR1
    File inputR2
    String outdirPath
    String? outPrefix = "flash"
    String? logPath = outdirPath + "/" + outPrefix + ".log"
    Int? minOverlap
    Int? maxOverlap
    Boolean? compress = true
    Int? threads = 1
    Int? memory = 4

    command {
        set -e -o pipefail
        mkdir -p ${outdirPath}
        ${preCommand}
        flash \
        ${"--threads=" + threads} \
        ${"--output-directory=" + outdirPath} \
        ${"--output-prefix=" + outPrefix} \
        ${true="--compress " false="" defined(compress)} \
        ${"--min-overlap=" + minOverlap} \
        ${"--max-overlap=" + maxOverlap} \
        ${inputR1} ${inputR2} 2>&1 | tee ${logPath}
    }

    output {
        File extendedFrags = outdirPath + "/" + outPrefix + ".extendedFrags.fastq.gz"
        File notCombined1 = outdirPath + "/" + outPrefix + ".notCombined_1.fastq.gz"
        File notCombined2 = outdirPath + "/" + outPrefix + ".notCombined_2.fastq.gz"
        File hist = outdirPath + "/" + outPrefix + ".hist"
        File histogram = outdirPath + "/" + outPrefix + ".histogram"
        String log = select_first([logPath])
    }

    runtime {
        cpu: select_first([threads])
        memory: select_first([memory])
    }

}