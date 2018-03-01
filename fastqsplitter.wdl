task FastqSplitter {
    File inputFastq
    String outputPath
    Int numberChunks
    File toolJar

    command {
    mkdir -p ${sep=' ' prefix(outputPath + "/chunk_", range(numberChunks))}
    ${if (numberChunks > 1) then ("java -jar " + toolJar + " -i " + inputFastq + write_lines(prefix("-o ", range(numberChunks))))
    else ("ln -s " + inputFastq + " " + outputPath + "/chunk_0/" + basename(inputFastq))}
    }

    output {
        Array[File] outputFastqFiles = glob(outputPath + "/chunk_*/" + basename(inputFastq))
    }
}