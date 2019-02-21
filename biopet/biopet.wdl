version 1.0

import "../common.wdl"

task BaseCounter {
    input {
        String? preCommand
        File? toolJar
        IndexedBamFile bam
        File refFlat
        String outputDir
        String prefix

        Int memory = 4
        Float memoryMultiplier = 3.5
    }

    String toolCommand = if defined(toolJar)
        then "java -Xmx" + memory + "G -jar " +toolJar
        else "biopet-basecounter -Xmx" + memory + "G"

    command {
        set -e -o pipefail
        mkdir -p ~{outputDir}
        ~{preCommand}
        ~{toolCommand} \
        -b ~{bam.file} \
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
        String dockerTag = "0.2--1"
    }

    String toolCommand = if defined(toolJar)
        then "java -Xmx" + memory + "G -jar " + toolJar
        else "biopet-extractadaptersfastqc -Xmx" + memory + "G"

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
        File adapterOutputFile = adapterOutputFilePath
        File contamsOutputFile = contamsOutputFilePath
        Array[String] adapterList = read_lines(adapterOutputFile)
        Array[String] contamsList = read_lines(contamsOutputFile)
    }

    runtime {
        memory: ceil(memory * memoryMultiplier)
        docker: "quay.io/biocontainers/biopet-extractadaptersfastqc:" + dockerTag
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
        FastqPair refFastq
        FastqPair inputFastq
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
        --in1 ~{inputFastq.R1} \
        --in2 ~{inputFastq.R2} \
        --ref1 ~{refFastq.R1} \
        --ref2 ~{refFastq.R2} \
        --out1 ~{out1path} \
        --out2 ~{out2path}
    }

    output {
        FastqPair out1 = object {
          R1: out1path,
          R1: out2path
        }
    }
    
    runtime {
        memory: ceil(memory * memoryMultiplier)
    }
}

task ReorderGlobbedScatters {
    input {
        Array[File]+ scatters
        # Should not be changed from the main pipeline. As it should not influence results.
        String dockerTag = "3.6"
    }

    command <<<
       set -e
       # Copy all the scatter files to the CWD so the output matches paths in
       # the cwd.
       for file in ~{sep=" " scatters}
          do cp $file .
       done
       python << CODE
       from os.path import basename
       scatters = ['~{sep="','" scatters}']
       splitext = [basename(x).split(".") for x in scatters]
       splitnum = [x.split("-") + [y] for x,y in splitext]
       ordered = sorted(splitnum, key=lambda x: int(x[1]))
       merged = ["{}-{}.{}".format(x[0],x[1],x[2]) for x in ordered]
       for x in merged:
           print(x)
       CODE
    >>>

    output {
        Array[File] reorderedScatters = read_lines(stdout())
    }

    runtime {
        docker: "python:" + dockerTag
        # 2 gigs of memory to be able to build the docker image in singularity
        memory: 2
    }
}

task ScatterRegions {
    input {
        Reference reference
        Int? scatterSize
        File? regions
        Boolean notSplitContigs = false
        Int memory = 4
        Float memoryMultiplier = 3.0
        String dockerTag = "0.2--0"
        String outputDirPath = "scatters"
    }

    command {
        set -e -o pipefail
        mkdir -p ~{outputDirPath}
        biopet-scatterregions -Xmx~{memory}G \
          -R ~{reference.fasta} \
          -o ~{outputDirPath} \
          ~{"-s " + scatterSize} \
          ~{"-L " + regions} \
          ~{true="--notSplitContigs" false="" notSplitContigs}
    }

    output {
        Array[File] scatters = glob(outputDirPath + "/scatter-*.bed")
    }

    runtime {
        docker: "quay.io/biocontainers/biopet-scatterregions:" + dockerTag
        memory: ceil(memory * memoryMultiplier)
    }
}

task ValidateAnnotation {
    input {
        String? preCommand
        File? toolJar
        File? refRefflat
        File? gtfFile
        Reference reference

        Int memory = 3
        Float memoryMultiplier = 3.0
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
        -R ~{reference.fasta}
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
        FastqPair inputFastq

        Int memory = 3
        Float memoryMultiplier = 3.0
    }

    String toolCommand = if defined(toolJar)
        then "java -Xmx" + memory + "G -jar " + toolJar
        else "biopet-validatefastq -Xmx" + memory + "G"

    command {
        set -e -o pipefail
        ~{preCommand}
        ~{toolCommand} \
        --fastq1 ~{inputFastq.R1} \
        ~{"--fastq2 " + inputFastq.R2}
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
        IndexedVcfFile vcf
        Reference reference

        Int memory = 3
        Float memoryMultiplier = 3.0
    }

    String toolCommand = if defined(toolJar)
        then "java -Xmx" + memory + "G -jar " + toolJar
        else "biopet-validatevcf -Xmx" + memory + "G"

    command {
        set -e -o pipefail
        ~{preCommand}
        ~{toolCommand} \
        -i ~{vcf.file} \
        -R ~{reference.fasta}
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
        IndexedVcfFile vcf
        Reference reference
        String outputDir
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
        Float memoryMultiplier = 2.5
        File? toolJar
        String? preCommand
    }

    String toolCommand = if defined(toolJar)
        then "java -Xmx" + memory + "G -jar " + toolJar
        else "biopet-vcfstats -Xmx" + memory + "G"

    command {
        set -e -o pipefail
        mkdir -p ~{outputDir}
        ~{preCommand}
        ~{toolCommand} \
        -I ~{vcf.file} \
        -R ~{reference.fasta} \
        -o ~{outputDir} \
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

    output {
        File? general = outputDir + "/general.tsv"
        File? genotype = outputDir + "/genotype.tsv"
        File? sampleDistributionAvailableAggregate = outputDir +
            "/sample_distributions/Available.aggregate.tsv"
        File? sampleDistributionAvailable = outputDir + "/sample_distributions/Available.tsv"
        File? sampleDistributionCalledAggregate = outputDir +
            "/sample_distributions/Called.aggregate.tsv"
        File? sampleDistributionCalled = outputDir + "/sample_distributions/Called.tsv"
        File? sampleDistributionFilteredAggregate = outputDir +
            "/sample_distributions/Filtered.aggregate.tsv"
        File? sampleDistributionFiltered = outputDir + "/sample_distributions/Filtered.tsv"
        File? sampleDistributionHetAggregate = outputDir + "/sample_distributions/Het.aggregate.tsv"
        File? sampleDistributionHetNoNRefAggregate = outputDir +
            "/sample_distributions/HetNonRef.aggregate.tsv"
        File? sampleDistributionHetNonRef = outputDir + "/sample_distributions/HetNonRef.tsv"
        File? sampleDistributionHet = outputDir + "/sample_distributions/Het.tsv"
        File? sampleDistributionHomAggregate = outputDir + "/sample_distributions/Hom.aggregate.tsv"
        File? sampleDistributionHomRefAggregate = outputDir +
            "/sample_distributions/HomRef.aggregate.tsv"
        File? sampleDistributionHomRef = outputDir + "/sample_distributions/HomRef.tsv"
        File? sampleDistributionHom = outputDir + "/sample_distributions/Hom.tsv"
        File? sampleDistributionHomVarAggregate = outputDir +
            "/sample_distributions/HomVar.aggregate.tsv"
        File? sampleDistributionHomVar = outputDir + "/sample_distributions/HomVar.tsv"
        File? sampleDistributionMixedAggregate = outputDir +
            "/sample_distributions/Mixed.aggregate.tsv"
        File? sampleDistributionMixed = outputDir + "/sample_distributions/Mixed.tsv"
        File? sampleDistributionNoCallAggregate = outputDir +
            "/sample_distributions/NoCall.aggregate.tsv"
        File? sampleDistributionNoCall = outputDir + "/sample_distributions/NoCall.tsv"
        File? sampleDistributionNonInformativeAggregate = outputDir +
            "/sample_distributions/NonInformative.aggregate.tsv"
        File? sampleDistributionNonInformative = outputDir +
            "/sample_distributions/NonInformative.tsv"
        File? sampleDistributionToalAggregate = outputDir +
            "/sample_distributions/Total.aggregate.tsv"
        File? sampleDistributionTotal = outputDir + "/sample_distributions/Total.tsv"
        File? sampleDistributionVariantAggregate = outputDir +
            "/sample_distributions/Variant.aggregate.tsv"
        File? sampleDistributionVariant = outputDir + "/sample_distributions/Variant.tsv"
        File? sampleCompareAlleleAbs = outputDir + "/sample_compare/allele.abs.tsv"
        File? sampleCompareAlleleNonRefAbs = outputDir + "/sample_compare/allele.non_ref.abs.tsv"
        File? sampleCompareAlleleRefAbs = outputDir + "/sample_compare/allele.ref.abs.tsv"
        File? sampleCompareAlleleRel = outputDir + "/sample_compare/allele.rel.tsv"
        File? sampleCompareGenotypeAbs = outputDir + "/sample_compare/genotype.abs.tsv"
        File? sampleCompareGenotypeNonRefAbs = outputDir +
            "/sample_compare/genotype.non_ref.abs.tsv"
        File? sampleCompareGenotypeRefAbs = outputDir + "/sample_compare/genotype.ref.abs.tsv"
        File? sampleCompareGenotypeRel = outputDir + "/sample_compare/genotype.rel.tsv"
    }

    runtime {
        cpu: localThreads
        memory: ceil(memory * memoryMultiplier)
    }
}
