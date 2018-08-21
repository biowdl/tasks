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

        Int memory = 4
        Float memoryMultiplier = 3.0
    }

    String toolCommand = if defined(toolJar)
        then "java -Xmx" + memory + "G -jar " +toolJar
        else "biopet-basecounter -Xmx" + memory + "G"

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
        memory: ceil(memory * memoryMultiplier)
    }
}

task ExtractAdaptersFastqc {
    input {
        File? toolJar
        String? preCommand
        File inputFile
        String outputDir
        String adapterOutputFilePath = outputDir + "/adapter.list"
        String contamsOutputFilePath = outputDir + "/contaminations.list"
        Boolean? skipContams
        File? knownContamFile
        File? knownAdapterFile
        Float? adapterCutoff
        Boolean? outputAsFasta

        Int memory = 4
        Float memoryMultiplier = 2.5
    }

    String toolCommand = if defined(toolJar)
        then "java -Xmx" + memory + "G -jar " +toolJar
        else "biopet-extractadaptersfastqc -Xmx" + memory + "G"

    command {
        set -e
        ~{preCommand}
        mkdir -p ~{outputDir}
        touch ~{adapterOutputFilePath}
        touch ~{contamsOutputFilePath}
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
        File adapterOutputFile = adapterOutputFilePath
        File contamsOutputFile = contamsOutputFilePath
        Array[String] adapterList = read_lines(adapterOutputFilePath)
        Array[String] contamsList = read_lines(contamsOutputFilePath)
    }

    runtime {
        memory: ceil(memory * memoryMultiplier)
    }
}

task FastqSplitter {
    input {
        String? preCommand
        File inputFastq
        Array[String]+ outputPaths
        File? toolJar

        Int memory = 4
        Float memoryMultiplier = 2.5
    }

    String toolCommand = if defined(toolJar)
        then "java -Xmx" + memory + "G -jar " +toolJar
        else "biopet-fastqsplitter -Xmx" + memory + "G"

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
        memory: ceil(memory * memoryMultiplier)
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

        Int memory = 4
        Float memoryMultiplier = 2.5
    }

    String toolCommand = if defined(toolJar)
        then "java -Xmx" + memory + "G -jar " +toolJar
        else "biopet-fastqsync -Xmx" + memory + "G"

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
        memory: ceil(memory * memoryMultiplier)
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

        Int memory = 4
        Float memoryMultiplier = 2.0
    }

    String toolCommand = if defined(toolJar)
        then "java -Xmx" + memory + "G -jar " +toolJar
        else "biopet-sampleconfig -Xmx" + memory + "G"

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
        memory: ceil(memory * memoryMultiplier)
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

        Int memory = 4
        Float memoryMultiplier = 3.0
    }

    String toolCommand = if defined(toolJar)
        then "java -Xmx" + memory + "G -jar " +toolJar
        else "biopet-scatterregions -Xmx" + memory + "G"

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
        memory: ceil(memory * memoryMultiplier)
    }
}

task Seqstat {
    input {
        String? preCommand
        File? toolJar
        File fastq
        String outputFile

        Int memory = 4
        Float memoryMultiplier = 2.0
    }

    String toolCommand = if defined(toolJar)
        then "java -Xmx" + memory + "G -jar " + toolJar
        else "biopet-seqstat -Xmx" + memory + "G"

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
        memory: ceil(memory * memoryMultiplier)
    }
}

task ValidateAnnotation {
    input {
        String? preCommand
        File? toolJar
        File? refRefflat
        File? gtfFile
        File refFasta
        File refFastaIndex
        File refDict

        Int memory = 4
        Float memoryMultiplier = 2.0
    }

    String toolCommand = if defined(toolJar)
        then "java -Xmx" + memory + "G -jar " + toolJar
        else "biopet-validateannotation -Xmx" + memory + "G"

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
        memory: ceil(memory * memoryMultiplier)
    }
}

task ValidateFastq {
    input {
        String? preCommand
        File? toolJar
        File fastq1
        File? fastq2

        Int memory = 4
        Float memoryMultiplier = 2.0
    }

    String toolCommand = if defined(toolJar)
        then "java -Xmx" + memory + "G -jar " + toolJar
        else "biopet-validatefastq -Xmx" + memory + "G"

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
        memory: ceil(memory * memoryMultiplier)
    }
}

task ValidateVcf {
    input {
        String? preCommand
        File? toolJar
        File vcfFile
        File vcfIndex
        File refFasta
        File refFastaIndex
        File refDict

        Int memory = 4
        Float memoryMultiplier = 2.0
    }

    String toolCommand = if defined(toolJar)
        then "java -Xmx" + memory + "G -jar " + toolJar
        else "biopet-validatevcf -Xmx" + memory + "G"

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
        memory: ceil(memory * memoryMultiplier)
    }
}

task VcfStats {
    input {
        File vcfFile
        File vcfIndex
        File refFasta
        File refFastaIndex
        File refDict
        String outDir
        File? intervals
        Array[String]+? infoTags
        Array[String]+? genotypeTags
        Int? sampleToSampleMinDepth
        Int? binSize
        Int? maxContigsInSingleJob
        Boolean writeBinStats = false
        Int localThreads = 1
        Boolean notWriteContigStats = false
        Boolean skipGeneral = false
        Boolean skipGenotype = false
        Boolean skipSampleDistributions = false
        Boolean skipSampleCompare = false
        String? sparkMaster
        Int? sparkExecutorMemory
        Array[String]+? sparkConfigValues

        Int memory = 4
        Float memoryMultiplier = 2.0
        File? toolJar
        String? preCommand
    }

    String toolCommand = if defined(toolJar)
        then "java -Xmx" + memory + "G -jar " + toolJar
        else "biopet-vcfstats -Xmx" + memory + "G"

    command {
        set -e -o pipefail
        ~{preCommand}
        ~{toolCommand} \
        -I ~{vcfFile} \
        -R ~{refFasta} \
        -o ~{outDir} \
        -t ~{localThreads} \
        ~{"--intervals " + intervals} \
        ~{true="--infoTag" false="" defined(infoTags)} ~{sep=" --infoTag " infoTags} \
        ~{true="--genotypeTag" false="" defined(genotypeTags)} ~{sep=" --genotypeTag "
            genotypeTags} \
        ~{"--sampleToSampleMinDepth " + sampleToSampleMinDepth} \
        ~{"--binSize " + binSize} \
        ~{"--maxContigsInSingleJob " + maxContigsInSingleJob} \
        ~{true="--writeBinStats" false="" writeBinStats} \
        ~{true="--notWriteContigStats" false="" notWriteContigStats} \
        ~{true="--skipGeneral" false="" skipGeneral} \
        ~{true="--skipGenotype" false="" skipGenotype} \
        ~{true="--skipSampleDistributions" false="" skipSampleDistributions} \
        ~{true="--skipSampleCompare" false="" skipSampleCompare} \
        ~{"--sparkMaster " + sparkMaster} \
        ~{"--sparkExecutorMemory " + sparkExecutorMemory} \
        ~{true="--sparkConfigValue" false="" defined(sparkConfigValues)} ~{
            sep=" --sparkConfigValue" sparkConfigValues}
    }

    runtime {
        cpu: localThreads
        memory: ceil(memory * memoryMultiplier)
    }
}
