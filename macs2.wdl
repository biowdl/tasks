version 1.0

task PeakCalling {
    input {
        Array[File]+ inputBams
        Array[File]+ inputBamsIndex
        Array[File]+? controlBams
        Array[File]+? controlBamsIndex
        String outDir
        String sampleName
        Boolean nomodel = false

        Int threads = 1
        Int memory = 8
    }

    command {
        set -e
        macs2 callpeak \
        --treatment ~{sep = ' ' inputBams} \
        ~{true="--control" false="" defined(controlBams)} ~{sep = ' ' controlBams} \
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
        docker: "quay.io/biocontainers/macs2:2.1.2--py27r351_0"
    }
}