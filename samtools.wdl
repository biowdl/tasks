task SamtoolsIndex {
    String bamFilePath

    command {
        samtools index ${bamFilePath}
    }

    output {
        File indexFile = bamFilePath + ".bai"
    }
}

task SamtoolsMerge {
    Array[File] bamFiles
    String outputBamPath

    command {
        samtools merge ${outputBamPath} ${bamFiles}
    }

    output {
        File bamFile = outputBamPath
    }
}