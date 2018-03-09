task FastqSplitter {
    File inputFastq
    String outputPath
    Int numberChunks
    String tool_jar
    Array[Int] chunks = range(numberChunks)

    command {
    mkdir -p ${sep=' ' prefix(outputPath + "/chunk_", chunks)}
    if [ ${numberChunks} -gt 1 ]; then
        SEP="/${basename(inputFastq)} -o "
        java -jar ${tool_jar} -I ${inputFastq} -o ${sep='$SEP' prefix(outputPath + "/chunk_", chunks)}/${basename(inputFastq)}
    else
        ln -sf ${inputFastq} ${outputPath}/chunk_0/${basename(inputFastq)}
    fi
    }

    output {
        Array[File] outputFastqFiles = glob(outputPath + "/chunk_*/" + basename(inputFastq))
    }
}

task ScatterRegions {
    File ref_fasta
    File ref_dict
    String outputDirPath
    String tool_jar
    Int? scatterSize
    File? regions

    command {
        mkdir -p ${outputDirPath}
        java -Xmx2G -jar ${tool_jar} \
          -R ${ref_fasta} \
          -o ${outputDirPath} \
          ${"-s " + scatterSize} \
          ${"-L " + regions}
    }

    output {
        Array[File] scatters = glob(outputDirPath + "/scatter-*.bed")
    }
}

task SampleConfig {
    String tool_jar
    Array[File]+ inputFiles
    String? sample
    String? library
    String? readgroup
    String? jsonOutputPath
    String? tsvOutputPath

    command {
        mkdir -p . $(dirname ${jsonOutputPath}) $(dirname ${tsvOutputPath})
        java -jar ${tool_jar} \
        -i ${sep="-i " inputFiles} \
        ${"--sample " + sample} \
        ${"--library " + library} \
        ${"--readgroup " + readgroup} \
        ${"--jsonOutput " + jsonOutputPath} \
        ${"--tsvOutput " + tsvOutputPath}
    }

    output {
        Array[String] keys = read_lines(stdout())
        File? jsonOutput = jsonOutputPath
        File? tsvOutput = tsvOutputPath
        Object values = if (defined(tsvOutput) && size(tsvOutput) > 0) then read_map(tsvOutput) else { "": "" }
    }
}

task DownloadSampleConfig {
    String? inputJar
    String? version = "0.1"

    command {
        ${if defined(inputJar) then "echo ${inputJar}" else "wget https://github.com/biopet/sampleconfig/releases/download/v${version}/SampleConfig-assembly-${version}.jar"}
    }

    output {
        File jar = if defined(inputJar) then select_first([inputJar]) else ("SampleConfig-assembly-" + version + ".jar")
    }
}
