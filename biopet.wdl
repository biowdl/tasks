version 1.0

task BaseCounter {
    input {
        String? preCommand
        File? toolJar
        File bam
        File bamIndex
        File refFlat
        String outputDir
        String prefix

        Float? memory
        Float? memoryMultiplier
    }

    Int mem = ceil(select_first([memory, 4.0]))
    String toolCommand = if defined(toolJar)
        then "java -Xmx" + mem + "G -jar " +toolJar
        else "biopet-basecounter -Xmx" + mem + "G"

    command {
        set -e -o pipefail
        mkdir -p ~{outputDir}
        ~{preCommand}
        ~{toolCommand} \
        -b ~{bam} \
        -r ~{refFlat} \
        -o ~{outputDir} \
        -p ~{prefix}
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

task ExtractAdaptersFastqc {
    input {
        File? toolJar
        String? preCommand
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
    }

    Int mem = ceil(select_first([memory, 4.0]))

    String toolCommand = if defined(toolJar)
        then "java -Xmx" + mem + "G -jar " +toolJar
        else "biopet-extractadaptersfastqc -Xmx" + mem + "G"

    command {
    set -e
    ~{preCommand}
    mkdir -p ~{outputDir}
    ~{toolCommand} \
    --inputFile ~{inputFile} \
    ~{"--adapterOutputFile " + adapterOutputFilePath } \
    ~{"--contamsOutputFile " + contamsOutputFilePath } \
    ~{"--knownContamFile " + knownContamFile} \
    ~{"--knownAdapterFile " + knownAdapterFile} \
    ~{"--adapterCutoff " + adapterCutoff} \
    ~{true="--skipContams" false="" skipContams} \
    ~{true="--outputAsFasta" false="" outputAsFasta}
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
    input {
        String? preCommand
        File inputFastq
        Array[String] outputPaths
        File? toolJar

        Float? memory
        Float? memoryMultiplier
    }

    Int mem = ceil(select_first([memory, 4.0]))
    String toolCommand = if defined(toolJar)
        then "java -Xmx" + mem + "G -jar " +toolJar
        else "biopet-fastqsplitter -Xmx" + mem + "G"

    command {
        set -e -o pipefail
        ~{preCommand}
        mkdir -p $(dirname ~{sep=') $(dirname ' outputPaths})
        if [ ~{length(outputPaths)} -gt 1 ]; then
            ~{toolCommand} \
            -I ~{inputFastq} \
            -o ~{sep=' -o ' outputPaths}
          else
            ln -sf ~{inputFastq} ~{outputPaths[0]}
          fi
    }

    output {
        Array[File] chunks = outputPaths
    }
    
    runtime {
        memory: ceil(mem * select_first([memoryMultiplier, 2.5]))
    }
}

task FastqSync {
    input {
        String? preCommand
        File ref1
        File ref2
        File in1
        File in2
        String out1path
        String out2path
        File? toolJar

        Float? memory
        Float? memoryMultiplier
    }

    Int mem = ceil(select_first([memory, 4.0]))

    String toolCommand = if defined(toolJar)
        then "java -Xmx" + mem + "G -jar " +toolJar
        else "biopet-fastqsync -Xmx" + mem + "G"

    command {
        set -e -o pipefail
        ~{preCommand}
        mkdir -p $(dirname ~{out1path}) $(dirname ~{out2path})
        ~{toolCommand} \
        --in1 ~{in1} \
        --in2 ~{in2} \
        --ref1 ~{ref1} \
        --ref2 ~{ref2} \
        --out1 ~{out1path} \
        --out2 ~{out2path}
    }

    output {
        File out1 = out1path
        File out2 = out2path
    }
    
    runtime {
        memory: ceil(mem * select_first([memoryMultiplier, 2.5]))
    }
}

task SampleConfig {
    input {
        File? toolJar
        String? preCommand
        Array[File]+ inputFiles
        String keyFilePath
        String? sample
        String? library
        String? readgroup
        String? jsonOutputPath
        String? tsvOutputPath
        Float? memory
        Float? memoryMultiplier
    }

    Int mem = ceil(select_first([memory, 4.0]))

    String toolCommand = if defined(toolJar)
        then "java -Xmx" + mem + "G -jar " +toolJar
        else "biopet-sampleconfig -Xmx" + mem + "G"

    command {
        set -e -o pipefail
        ~{preCommand}
        mkdir -p . ~{"$(dirname " + jsonOutputPath + ")"} ~{"$(dirname " + tsvOutputPath + ")"}
        ~{toolCommand} \
        -i ~{sep="-i " inputFiles} \
        ~{"--sample " + sample} \
        ~{"--library " + library} \
        ~{"--readgroup " + readgroup} \
        ~{"--jsonOutput " + jsonOutputPath} \
        ~{"--tsvOutput " + tsvOutputPath} \
        > ~{keyFilePath}
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
    input {
        String? preCommand
        File refFasta
        File refDict
        String outputDirPath
        File? toolJar
        Int? scatterSize
        File? regions

        Float? memory
        Float? memoryMultiplier
    }
    Int mem = ceil(select_first([memory, 4.0]))

    String toolCommand = if defined(toolJar)
        then "java -Xmx" + mem + "G -jar " +toolJar
        else "biopet-scatterregions -Xmx" + mem + "G"

    command {
        set -e -o pipefail
        ~{preCommand}
        mkdir -p ~{outputDirPath}
        ~{toolCommand} \
          -R ~{refFasta} \
          -o ~{outputDirPath} \
          ~{"-s " + scatterSize} \
          ~{"-L " + regions}
    }

    output {
        Array[File] scatters = glob(outputDirPath + "/scatter-*.bed")
    }

    runtime {
        memory: ceil(mem * select_first([memoryMultiplier, 3.0]))
    }
}

task Seqstat {
    input {
        String? preCommand
        File? toolJar
        File fastq
        String outputFile
        Float? memory
        Float? memoryMultiplier
    }

    Int mem = ceil(select_first([memory, 4.0]))

    String toolCommand = if defined(toolJar)
        then "java -Xmx" + mem + "G -jar " + toolJar
        else "biopet-seqstat -Xmx" + mem + "G"

    command {
        set -e -o pipefail
        ~{preCommand}
        mkdir -p $(dirname ~{outputFile})
        ~{toolCommand} \
        --fastq ~{fastq} \
        --output ~{outputFile}
    }

    output {
        File json = outputFile
    }

    runtime {
        memory: ceil(mem * select_first([memoryMultiplier, 2.0]))
    }
}

task ValidateAnnotation {
    input {
        String? preCommand
        File? toolJar
        File? refRefflat
        File? gtfFile
        File refFasta
        Float? memory
        Float? memoryMultiplier
    }

    Int mem = ceil(select_first([memory, 4.0]))

    String toolCommand = if defined(toolJar)
        then "java -Xmx" + mem + "G -jar " + toolJar
        else "biopet-validateannotation -Xmx" + mem + "G"

    command {
        set -e -o pipefail
        ~{preCommand}
        ~{toolCommand} \
        ~{"-r " + refRefflat} \
        ~{"-g " + gtfFile} \
        -R ~{refFasta}
    }

    output {
        File stderr = stderr()
    }

    runtime {
        memory: ceil(mem * select_first([memoryMultiplier, 2.0]))
    }
}

task ValidateFastq {
    input {
        String? preCommand
        File? toolJar
        File fastq1
        File? fastq2

        Float? memory
        Float? memoryMultiplier
    }

    Int mem = ceil(select_first([memory, 4.0]))

    String toolCommand = if defined(toolJar)
        then "java -Xmx" + mem + "G -jar " + toolJar
        else "biopet-validatefastq -Xmx" + mem + "G"

    command {
        set -e -o pipefail
        ~{preCommand}
        ~{toolCommand} \
        --fastq1 ~{fastq1} \
        ~{"--fastq2 " + fastq2}
    }

    output {
        File stderr = stderr()
    }

    runtime {
        memory: ceil(mem * select_first([memoryMultiplier, 2.0]))
    }
}

task ValidateVcf {
    input {
        String? preCommand
        File? toolJar
        File vcfFile
        File refFasta

        Float? memory
        Float? memoryMultiplier
    }

    Int mem = ceil(select_first([memory, 4.0]))

    String toolCommand = if defined(toolJar)
        then "java -Xmx" + mem + "G -jar " + toolJar
        else "biopet-validatevcf -Xmx" + mem + "G"

    command {
        set -e -o pipefail
        ~{preCommand}
        ~{toolCommand} \
        -i ~{vcfFile} \
        -R ~{refFasta}
    }

    output {
        File stderr = stderr()
    }

    runtime {
        memory: ceil(mem * select_first([memoryMultiplier, 2.0]))
    }
}
