#workflow runFlash {
#    String preCommand
#    String outDirPath
#    File R1
#    File R2
#
#    call flash {
#        input:
#            preCommand = preCommand,
#            outDirPath = outDirPath,
#            inputR1 = R1,
#            inputR2 = R2,
#    }
#}
task flash {
    String? preCommand
    File inputR1
    File inputR2
    String outDirPath
    String? outPrefix = "flash"
    String? logPath = outDirPath + "/" + outPrefix + ".log"
    Int? minOverlap
    Int? maxOverlap
    Boolean? compress = true
    Int? threads = 1
    Int? memory = 4

    command {
        set -e -o pipefail
        mkdir -p ${outDirPath}
        ${preCommand}
        flash \
        ${"--threads=" + threads} \
        ${"--output-directory=" + outDirPath} \
        ${"--output-prefix=" + outPrefix} \
        ${true="--compress " false="" defined(compress)} \
        ${"--min-overlap=" + minOverlap} \
        ${"--max-overlap=" + maxOverlap} \
        ${inputR1} ${inputR2} 2>&1 | tee ${logPath}
    }

    output {
        File extendedFrags = outDirPath + "/" + outPrefix + ".extendedFrags.fastq.gz"
        File notCombined1 = outDirPath + "/" + outPrefix + ".notCombined_1.fastq.gz"
        File notCombined2 = outDirPath + "/" + outPrefix + ".notCombined_2.fastq.gz"
        File hist = outDirPath + "/" + outPrefix + ".hist"
        File histogram = outDirPath + "/" + outPrefix + ".histogram"
        String log = select_first([logPath])
    }

    runtime {
        cpu: select_first([threads])
        memory: select_first([memory])
    }

}