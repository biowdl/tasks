task PeakCalling {
    String? preCommand
    Array[File] bamFiles
    String outDir
    String sampleName
    Int? threads
    Int? memory

    command {
        set -e -o pipefail
        ${preCommand}
        macs2 callpeak \
        --treatment ${sep = ' ' bamFiles} \
        --outdir ${outDir} \
        --name ${sampleName}
    }

    output {
        File peakFile = outDir + "/" + sampleName + "/macs2/" + sampleName + "_peaks.narrowPeak"
    }

    runtime {
        cpu: select_first([threads,1])
        memory: select_first([memory,8])
    }
}