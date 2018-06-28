task peakCalling {
    String? preCommand
    File bamFile
    String outDir
    String sampleName
    Int? threads
    Int? memory


    command {
        set -e -o pipefail
        ${preCommand}
        macs2 callpeaks --treatment ${bamFile} --outdir ${outDir} --name ${sampleName}
    }

    output {
        File peakFile = outDir + "/" + sampleName + "/macs2/" + sampleName + "_peaks.narrowPeakd"
    }

    runtime {
        cpu: select_first([threads,1])
        memory: select_first([memory,8])
    }
}