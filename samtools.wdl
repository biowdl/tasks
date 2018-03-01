task Index {
    String bamFilePath

    command {
        samtools index ${bamFilePath}
    }

    output {
        File indexFile = bamFilePath + ".bai"
    }
}

task Merge {
    Array[File]+ bamFiles
    String outputBamPath

    command {
        samtools merge ${outputBamPath} ${sep=' ' bamFiles}
    }

    output {
        File outputBam = outputBamPath
    }
}