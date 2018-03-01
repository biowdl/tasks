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

task Markdup {
    File inputBam
    String outputBamPath

    command {
        samtools markdup ${inputBam} ${outputBamPath}
    }

    output {
        File outputBam = outputBamPath
    }
}

task Flagstat {
    File inputBam
    String outputPath

    command {
        samtools flagstat ${inputBam} > ${outputPath}
    }

    output {
        File flagstat = outputPath
    }
}
