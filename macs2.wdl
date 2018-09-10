version 1.0

import "common.wdl"

task PeakCalling {
    input {
        String? preCommand
        Array[IndexedBamFile] bamFiles
        String outDir
        String sampleName
        Int threads = 1
        Int memory = 8
        Boolean nomodel = false
    }

    command {
        set -e -o pipefail
        ~{preCommand}
        macs2 callpeak \
        --treatment ~{sep = ' ' bamFiles.file} \
        --outdir ~{outDir} \
        --name ~{sampleName} \
        ~{true='--nomodel' false='' nomodel}
    }

    output {
        File peakFile = outDir + "/" + sampleName + "_peaks.narrowPeak"
    }

    runtime {
        cpu: threads
        memory: memory
    }
}