version 1.0

task PeakCalling {
    input {
        String? preCommand
        Array[File] bamFiles
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
        --treatment ~{sep = ' ' bamFiles} \
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