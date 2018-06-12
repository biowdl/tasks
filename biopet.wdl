# PLEASE ADD TASKS IN ALPHABETIC ORDER.
# This makes searching a lot easier.

task BaseCounter {
    String? preCommand
    File toolJar
    File bam
    File bamIndex
    File refFlat
    String outputDir
    String prefix

    Float? memory
    Float? memoryMultiplier

    Int mem = ceil(select_first([memory, 12.0]))
    command {
        set -e -o pipefail
        mkdir -p ${outputDir}
        ${preCommand}
        java -Xmx${mem}G -jar ${toolJar} \
        -b ${bam} \
        -r ${refFlat} \
        -o ${outputDir} \
        -p ${prefix}
    }

    output {
        File exonAntisense = outputDir + "/" + prefix + ".base.exon.antisense.counts"
        File exon = outputDir + "/" + prefix + ".base.exon.counts"
        File exonMergeAntisense = outputDir + "/" + prefix + ".base.exon.merge.antisense.counts"
        File exonMerge = outputDir + "/" + prefix + ".base.exon.merge.counts"
        File exonMergeSense = outputDir + "/" + prefix + ".base.exon.merge.sense.counts"
        File exonSense = outputDir + "/" + prefix + ".base.exon.sense.counts"
        File geneAntisense = outputDir + "/" + prefix + ".base.gene.antisense.counts"
        File gene = outputDir + "/" + prefix + ".base.gene.counts"
        File geneExonicAntisense = outputDir + "/" + prefix + ".base.gene.exonic.antisense.counts"
        File geneExonic = outputDir + "/" + prefix + ".base.gene.exonic.counts"
        File geneExonicSense = outputDir + "/" + prefix + ".base.gene.exonic.sense.counts"
        File geneIntronicAntisense = outputDir + "/" + prefix + ".base.gene.intronic.antisense.counts"
        File geneIntronic = outputDir + "/" + prefix + ".base.gene.intronic.counts"
        File geneIntronicSense = outputDir + "/" + prefix + ".base.gene.intronic.sense.counts"
        File geneSense = outputDir + "/" + prefix + ".base.gene.sense.counts"
        File intronAntisense = outputDir + "/" + prefix + ".base.intron.antisense.counts"
        File intron = outputDir + "/" + prefix + ".base.intron.counts"
        File intronMergeAntisense = outputDir + "/" + prefix + ".base.intron.merge.antisense.counts"
        File intronMerge = outputDir + "/" + prefix + ".base.intron.merge.counts"
        File intronMergeSense = outputDir + "/" + prefix + ".base.intron.merge.sense.counts"
        File intronSense = outputDir + "/" + prefix + ".base.intron.sense.counts"
        File metaExonsNonStranded = outputDir + "/" + prefix + ".base.metaexons.non_stranded.counts"
        File metaExonsStrandedAntisense = outputDir + "/" + prefix + ".base.metaexons.stranded.antisense.counts"
        File metaExonsStranded = outputDir + "/" + prefix + ".base.metaexons.stranded.counts"
        File metaExonsStrandedSense = outputDir + "/" + prefix + ".base.metaexons.stranded.sense.counts"
        File transcriptAntisense = outputDir + "/" + prefix + ".base.transcript.antisense.counts"
        File transcript = outputDir + "/" + prefix + ".base.transcript.counts"
        File transcriptExonicAntisense = outputDir + "/" + prefix + ".base.transcript.exonic.antisense.counts"
        File transcriptExonic = outputDir + "/" + prefix + ".base.transcript.exonic.counts"
        File transcriptExonicSense = outputDir + "/" + prefix + ".base.transcript.exonic.sense.counts"
        File transcriptIntronicAntisense = outputDir + "/" + prefix + ".base.transcript.intronic.antisense.counts"
        File transcriptIntronic = outputDir + "/" + prefix + ".base.transcript.intronic.counts"
        File transcriptIntronicSense = outputDir + "/" + prefix + ".base.transcript.intronic.sense.counts"
        File transcriptSense = outputDir + "/" + prefix + ".base.transcript.sense.counts"
    }

    runtime {
        memory: ceil(mem * select_first([memoryMultiplier, 3.0]))
    }
}

task extractAdaptersFastqc {
    File? toolJar
    File inputFile
    String outputDir
    String? adapterOutputFilePath = outputDir + "/adapter.list"
    String? contamsOutputFilePath = outputDir + "/contaminations.list"
    Boolean? skipContams
    File? knownContamFile
    File? knownAdapterFile
    Float? adapterCutoff
    Boolean? outputAsFasta

    Float? memory
    Float? memoryMultiplier
    Int mem = ceil(select_first([memory, 4.0]))

    String toolCommand = if defined(toolJar)
    then "java -Xmx" + mem + "G -jar " +toolJar
    else "biopet-extractadaptersfastqc -Xmx" + mem + "G"

    command {
    set -e
    mkdir -p ${outputDir}
    ${toolCommand} \
    --inputFile ${inputFile} \
    ${"--adapterOutputFile " + adapterOutputFilePath } \
    ${"--contamsOutputFile " + contamsOutputFilePath } \
    ${"--knownContamFile " + knownContamFile} \
    ${"--knownAdapterFile " + knownAdapterFile} \
    ${"--adapterCutoff " + adapterCutoff} \
    ${true="--skipContams" false="" skipContams} \
    ${true="--outputAsFasta" false="" outputAsFasta}
    }

    output {
        File adapterOutputFile = select_first([adapterOutputFilePath])
        File contamsOutputFile = select_first([contamsOutputFilePath])
        Array[String] adapterList = read_lines(select_first([adapterOutputFilePath]))
        Array[String] contamsList = read_lines(select_first([contamsOutputFilePath]))
    }

    runtime {
        memory: ceil(mem * select_first([memoryMultiplier, 2.5]))
    }
}


task FastqSplitter {
    String? preCommand
    File inputFastq
    String outputPath
    Int numberChunks
    String tool_jar
    Array[Int] chunks = range(numberChunks)

    command {
        set -e -o pipefail
        ${preCommand}
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

task FastqSync {
    String? preCommand
    File ref1
    File ref2
    File in1
    File in2
    String out1path
    String out2path
    File tool_jar
    command {
        set -e -o pipefail
        ${preCommand}
        mkdir -p $(dirname ${out1path}) $(dirname ${out2path})
        java -jar ${tool_jar} \
        --in1 ${in1} \
        --in2 ${in2} \
        --ref1 ${ref1} \
        --ref2 ${ref2} \
        --out1 ${out1path} \
        --out2 ${out2path}
    }
    output {
        File out1 = out1path
        File out2 = out2path
    }
}

task SampleConfig {
    String? preCommand
    String tool_jar
    Array[File]+ inputFiles
    String keyFilePath
    String? sample
    String? library
    String? readgroup
    String? jsonOutputPath
    String? tsvOutputPath

    Float? memory
    Float? memoryMultiplier

    Int mem = ceil(select_first([memory, 4.0]))
    command {
        set -e -o pipefail
        ${preCommand}
        mkdir -p . ${"$(dirname " + jsonOutputPath + ")"} ${"$(dirname " + tsvOutputPath + ")"}
        java -Xmx${mem}G -jar ${tool_jar} \
        -i ${sep="-i " inputFiles} \
        ${"--sample " + sample} \
        ${"--library " + library} \
        ${"--readgroup " + readgroup} \
        ${"--jsonOutput " + jsonOutputPath} \
        ${"--tsvOutput " + tsvOutputPath} \
        > ${keyFilePath}
    }

    output {
        File keysFile = keyFilePath
        File? jsonOutput = jsonOutputPath
        File? tsvOutput = tsvOutputPath
    }

    runtime {
        memory: ceil(mem * select_first([memoryMultiplier, 2.0]))
    }
}

task ScatterRegions {
    String? preCommand
    File ref_fasta
    File ref_dict
    String outputDirPath
    String tool_jar
    Int? scatterSize
    File? regions

    Float? memory
    Float? memoryMultiplier

    Int mem = ceil(select_first([memory, 4.0]))
    command {
        set -e -o pipefail
        ${preCommand}
        mkdir -p ${outputDirPath}
        java -Xmx${mem}G -jar ${tool_jar} \
          -R ${ref_fasta} \
          -o ${outputDirPath} \
          ${"-s " + scatterSize} \
          ${"-L " + regions}
    }

    output {
        Array[File] scatters = glob(outputDirPath + "/scatter-*.bed")
    }

    runtime {
        memory: ceil(mem * select_first([memoryMultiplier, 3.0]))
    }
}

