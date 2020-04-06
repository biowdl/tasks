version 1.0

# Copyright (c) 2017 Leiden University Medical Center
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import "../common.wdl"

task BaseCounter {
    input {
        String? preCommand
        File? toolJar
        IndexedBamFile bam
        File refFlat
        String outputDir
        String prefix

        String memory = "14G"
        String javaXmx = "4G"
    }

    String toolCommand = if defined(toolJar)
        then "java -Xmx~{javaXmx} -jar " + toolJar
        else "biopet-basecounter -Xmx~{javaXmx}"

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
        memory: memory
    }
}

task ExtractAdaptersFastqc {
    input {
        File inputFile
        String outputDir
        String adapterOutputFilePath = outputDir + "/adapter.list"
        String contamsOutputFilePath = outputDir + "/contaminations.list"
        Boolean? skipContams
        File? knownContamFile
        File? knownAdapterFile
        Float? adapterCutoff
        Boolean? outputAsFasta

        String memory = "40G" # This is ridiculous, but needed due to vmem monitoring on SGE.
        String javaXmx = "8G"
        String dockerImage = "quay.io/biocontainers/biopet-extractadaptersfastqc:0.2--1"
    }

    command {
        set -e
        mkdir -p ~{outputDir}
        biopet-extractadaptersfastqc -Xmx~{javaXmx} \
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
        memory: memory
        docker: dockerImage
    }
}

task FastqSplitter {
    input {
        String? preCommand
        File inputFastq
        Array[String]+ outputPaths
        File? toolJar

        String memory = "12G"
        String javaXmx = "4G"
        String dockerImage = "quay.io/biocontainers/biopet-fastqsplitter:0.1--2"
    }

    command {
        set -e
        mkdir -p $(dirname ~{sep=') $(dirname ' outputPaths})
        biopet-fastqsplitter -Xmx~{javaXmx} \
        -I ~{inputFastq} \
        -o ~{sep=' -o ' outputPaths}
    }

    output {
        Array[File] chunks = outputPaths
    }
    
    runtime {
        memory: memory
        docker: dockerImage
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

        String memory = "10G"
        String javaXmx = "4G"
    }

    String toolCommand = if defined(toolJar)
        then "java -Xmx~{javaXmx} -jar " + toolJar
        else "biopet-fastqsync -Xmx~{javaXmx}"

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
          R2: out2path
        }
    }
    
    runtime {
        memory: memory
    }
}

task ReorderGlobbedScatters {
    input {
        Array[File]+ scatters

        # Should not be changed from the main pipeline. As it should not influence results.
        # The 3.7-slim container is 143 mb on the filesystem. 3.7 is 927 mb.
        # The slim container is sufficient for this small task.
        String dockerImage = "python:3.7-slim"
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
        docker: dockerImage
        # 4 gigs of memory to be able to build the docker image in singularity
        memory: "4G"
    }

    parameter_meta {
        scatters: {description: "The files which should be ordered.", category: "required"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}

task ScatterRegions {
    input {
        File referenceFasta
        File referenceFastaDict
        Int? scatterSize
        File? regions
        Boolean notSplitContigs = false
        File? bamFile
        File? bamIndex

        String memory = "1G"
        String javaXmx = "500M"
        Int timeMinutes = 1
        String dockerImage = "quay.io/biocontainers/biopet-scatterregions:0.2--0"
    }

    # OutDirPath must be defined here because the glob process relies on
    # linking. This path must be in the containers filesystem, otherwise the
    # linking does not work.
    String outputDirPath = "scatters"

    command <<<
        set -e -o pipefail
        mkdir -p ~{outputDirPath}
        biopet-scatterregions -Xmx~{javaXmx} -XX:ParallelGCThreads=1 \
          -R ~{referenceFasta} \
          -o ~{outputDirPath} \
          ~{"-s " + scatterSize} \
          ~{"-L " + regions} \
          ~{"--bamFile " + bamFile} \
          ~{true="--notSplitContigs" false="" notSplitContigs}

        # Glob messes with order of scatters (10 comes before 1), which causes
        # problems at gatherGvcfs
        # Therefore we reorder the scatters with python.
        python << CODE
        import os
        scatters = os.listdir("~{outputDirPath}")
        splitext = [ x.split(".") for x in scatters]
        splitnum = [x.split("-") + [y] for x,y in splitext]
        ordered = sorted(splitnum, key=lambda x: int(x[1]))
        merged = ["~{outputDirPath}/{}-{}.{}".format(x[0],x[1],x[2]) for x in ordered]
        for x in merged:
          print(x)
        CODE
    >>>

    output {
        Array[File] scatters =  read_lines(stdout())
    }

    runtime {
        docker: dockerImage
        time_minutes: timeMinutes
        memory: memory
    }

    parameter_meta {
        referenceFasta: {description: "The reference fasta file.", category: "required"}
        referenceFastaDict: {description: "The sequence dictionary associated with the reference fasta file.",
                             category: "required"}
        scatterSize: {description: "Equivalent to biopet scatterregions' `-s` option.", category: "common"}
        regions: {description: "The regions to be scattered.", category: "advanced"}
        notSplitContigs: {description: "Equivalent to biopet scatterregions' `--notSplitContigs` flag.",
                          category: "advanced"}
        bamFile: {description: "Equivalent to biopet scatterregions' `--bamfile` option.",
                  category: "advanced"}
        bamIndex: {description: "The index for the bamfile given through bamFile.", category: "advanced"}

        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.",
                  category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}

task ValidateAnnotation {
    input {
        File? refRefflat
        File? gtfFile
        Reference reference

        String memory = "9G"
        String javaXmx = "3G"
        String dockerImage = "quay.io/biocontainers/biopet-validateannotation:0.1--0"
    }

    command {
        biopet-validateannotation -Xmx~{javaXmx} \
        ~{"-r " + refRefflat} \
        ~{"-g " + gtfFile} \
        -R ~{reference.fasta}
    }

    output {
        File stderr = stderr()
    }

    runtime {
        memory: memory
        docker: dockerImage
    }
}

task ValidateFastq {
    input {
        File read1
        File? read2
        String memory = "9G"
        String javaXmx = "3G"
        String dockerImage = "quay.io/biocontainers/biopet-validatefastq:0.1.1--1"
    }

    command {
        biopet-validatefastq -Xmx~{javaXmx} \
        --fastq1 ~{read1} \
        ~{"--fastq2 " + read2}
    }

    output {
        File stderr = stderr()
    }

    runtime {
        memory: memory
        docker: dockerImage
    }
}

task ValidateVcf {
    input {
        IndexedVcfFile vcf
        Reference reference
        String memory = "9G"
        String javaXmx = "3G"
        String dockerImage = "quay.io/biocontainers/biopet-validatevcf:0.1--0"
    }

    command {
        biopet-validatevcf -Xmx~{javaXmx} \
        -i ~{vcf.file} \
        -R ~{reference.fasta}
    }

    output {
        File stderr = stderr()
    }

    runtime {
        memory: memory
        docker: dockerImage
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

        String dockerImage = "quay.io/biocontainers/biopet-vcfstats:1.2--0"
        String memory = "12G"
        String javaXmx = "4G"
    }

    command {
        set -e
        mkdir -p ~{outputDir}
        biopet-vcfstats -Xmx~{javaXmx} \
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
        # A glob is easier, but duplicates all the outputs
        Array[File] allStats = select_all([
            general,
            genotype,
            sampleDistributionAvailableAggregate,
            sampleDistributionAvailable,
            sampleDistributionCalledAggregate,
            sampleDistributionCalled,
            sampleDistributionFilteredAggregate,
            sampleDistributionFiltered,
            sampleDistributionHetAggregate,
            sampleDistributionHetNoNRefAggregate,
            sampleDistributionHetNonRef,
            sampleDistributionHet,
            sampleDistributionHomAggregate,
            sampleDistributionHomRefAggregate,
            sampleDistributionHomRef,
            sampleDistributionHom,
            sampleDistributionHomVarAggregate,
            sampleDistributionHomVar,
            sampleDistributionMixedAggregate,
            sampleDistributionMixed,
            sampleDistributionNoCallAggregate,
            sampleDistributionNoCall,
            sampleDistributionNonInformativeAggregate,
            sampleDistributionNonInformative,
            sampleDistributionToalAggregate,
            sampleDistributionTotal,
            sampleDistributionVariantAggregate,
            sampleDistributionVariant,
            sampleCompareAlleleAbs,
            sampleCompareAlleleNonRefAbs,
            sampleCompareAlleleRefAbs,
            sampleCompareAlleleRel,
            sampleCompareGenotypeAbs,
            sampleCompareGenotypeNonRefAbs,
            sampleCompareGenotypeRefAbs,
            sampleCompareGenotypeRel
        ])
    }

    runtime {
        cpu: localThreads
        memory: memory
        docker: dockerImage
    }
}
