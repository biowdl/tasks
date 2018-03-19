task MergeCounts {
    Array[File] inputFiles
    String outputFile
    String idVar
    String measurementVar
    File script

    command {
        Rscript ${script} \
        ${idVar} \
        ${measurementVar} \
        ${sep=" " inputFiles} \
        > ${outputFile}
    }

    output {
        File mergedCounts = outputFile
    }
}