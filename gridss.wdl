version 1.0

# Copyright (c) 2020 Leiden University Medical Center
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

import "bwa.wdl" as bwa

task AnnotateInsertedSequence {
    input {
        File inputVcf
        String outputPath = "gridss.annotated.vcf.gz"
        File viralReference
        File viralReferenceFai
        File viralReferenceDict
        File viralReferenceImg

        Int threads = 8
        String javaXmx = "8G"
        String memory = "9G"
        String dockerImage = "quay.io/biowdl/gridss:2.12.2"
        Int timeMinutes = 120
    }

    command {
        set -e
        _JAVA_OPTIONS="$_JAVA_OPTIONS -Xmx~{javaXmx}"
        AnnotateInsertedSequence \
        REFERENCE_SEQUENCE=~{viralReference} \
        INPUT=~{inputVcf} \
        OUTPUT=~{outputPath} \
        ALIGNMENT=APPEND \
        WORKING_DIR='.' \
        WORKER_THREADS=~{threads}
    }

    output {
        File outputVcf = outputPath
        File outputVcfIndex = outputPath + ".tbi"
    }

    runtime {
        cpu: threads
        memory: memory
        time_minutes: timeMinutes # !UnknownRuntimeKey
        docker: dockerImage
    }

    parameter_meta {
        inputVcf: {description: "The input VCF file.", category: "required"}
        outputPath: {description: "The path the output will be written to.", category: "common"}
        viralReference: {description: "A fasta file with viral sequences.", category: "required"}
        viralReferenceFai: {description: "The index for the viral reference fasta.", category: "required"}
        viralReferenceDict: {description: "The dict file for the viral reference.", category: "required"}
        viralReferenceImg: {description: "The BWA index image (generated with GATK BwaMemIndexImageCreator) of the viral reference.", category: "required"}

        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.",
                  category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}

task AnnotateSvTypes {
    input {
        File gridssVcf
        File gridssVcfIndex
        String outputPath = "./gridss.svtyped.vcf.bgz"

        String memory = "32G"
        String dockerImage = "quay.io/biocontainers/bioconductor-structuralvariantannotation:1.10.0--r41hdfd78af_0"
        Int timeMinutes = 240
    }

    String effectiveOutputPath = sub(outputPath, "\\.bgz", "")
    String index = if effectiveOutputPath != outputPath then "T" else "F"


    # Based on https://github.com/PapenfussLab/gridss/issues/74
    command <<<
        set -e
        mkdir -p "$(dirname ~{outputPath})"
        R --vanilla << "EOF"
        library(VariantAnnotation)
        library(StructuralVariantAnnotation)

        vcf_path <- "~{gridssVcf}"
        out_path <- "~{effectiveOutputPath}"

        # Simple SV type classifier
        simpleEventType <- function(gr) {
          return(ifelse(seqnames(gr) != seqnames(partner(gr)), "BND", # inter-chromosomosal
                  ifelse(gr$insLen >= abs(gr$svLen) * 0.7, "INS",
                   ifelse(strand(gr) == strand(partner(gr)), "INV",
                    ifelse(xor(start(gr) < start(partner(gr)), strand(gr) == "-"), "DEL",
                     "DUP")))))
        }

        header <- scanVcfHeader(vcf_path)
        vcf <- readVcf(vcf_path, seqinfo(header))
        gr <- breakpointRanges(vcf)
        svtype <- simpleEventType(gr)
        info(vcf[gr$sourceId])$SVTYPE <- svtype
        # GRIDSS doesn't supply a GT, so we estimate GT based on AF (assuming CN of 2, might be inaccurate)
        geno(vcf)$GT <- ifelse(geno(vcf)$AF > 0.75, "1/1", ifelse(geno(vcf)$AF < 0.25, "0/0", "0/1"))
        writeVcf(vcf, out_path, index=~{index})
        EOF
    >>>

    output {
        File vcf = outputPath
        File? vcfIndex = outputPath + ".tbi"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes # !UnknownRuntimeKey
        docker: dockerImage
    }

    parameter_meta {
        gridssVcf: {description: "The VCF produced by GRIDSS.", category: "required"}
        gridssVcfIndex: {description: "The index for the VCF produced by GRIDSS.", category: "required"}
        outputPath: {description: "The path the output should be written to.",  category: "common"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}

task GRIDSS {
    input {
        File tumorBam
        File tumorBai
        String tumorLabel
        BwaIndex reference
        String outputPrefix = "gridss"

        File? normalBam
        File? normalBai
        String? normalLabel
        File? blacklistBed
        File? gridssProperties

        Int jvmHeapSizeGb = 64
        Int nonJvmMemoryGb = 10
        Int threads = 12
        Int timeMinutes = ceil(7200 / threads) + 1800
        String dockerImage = "quay.io/biowdl/gridss:2.12.2"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPrefix})"
        gridss \
        -w . \
        --reference ~{reference.fastaFile} \
        --output ~{outputPrefix}.vcf.gz \
        --assembly ~{outputPrefix}_assembly.bam \
        ~{"-c " + gridssProperties} \
        ~{"-t " + threads} \
        ~{"--jvmheap " + jvmHeapSizeGb + "G"} \
        --labels ~{normalLabel}~{true="," false="" defined(normalLabel)}~{tumorLabel} \
        ~{"--blacklist " + blacklistBed} \
        ~{normalBam} \
        ~{tumorBam}
        samtools index ~{outputPrefix}_assembly.bam ~{outputPrefix}_assembly.bai

        # For some reason the VCF index is sometimes missing
        if [ ! -e ~{outputPrefix}.vcf.gz.tbi ]
          then
            tabix ~{outputPrefix}.vcf.gz
        fi
    }

    output {
        File vcf = outputPrefix + ".vcf.gz"
        File vcfIndex = outputPrefix + ".vcf.gz.tbi"
        File assembly = outputPrefix + "_assembly.bam"
        File assemblyIndex = outputPrefix + "_assembly.bai"
    }

    runtime {
        cpu: threads
        memory: "~{jvmHeapSizeGb + nonJvmMemoryGb}G"
        time_minutes: timeMinutes # !UnknownRuntimeKey
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        tumorBam: {description: "The input BAM file. This should be the tumor/case sample in case of a paired analysis.", category: "required"}
        tumorBai: {description: "The index for tumorBam.", category: "required"}
        tumorLabel: {description: "The name of the (tumor) sample.", category: "required"}
        reference: {description: "A BWA index, this should also include the fasta index file (.fai).", category: "required"}
        outputPrefix: {description: "The prefix for the output files. This may include parent directories.", category: "common"}
        normalBam: {description: "The BAM file for the normal/control sample.", category: "advanced"}
        normalBai: {description: "The index for normalBam.", category: "advanced"}
        normalLabel: {description: "The name of the normal sample.", category: "advanced"}
        blacklistBed: {description: "A bed file with blaclisted regins.", category: "advanced"}
        gridssProperties: {description: "A properties file for gridss.", category: "advanced"}

        threads: {description: "The number of the threads to use.", category: "advanced"}
        jvmHeapSizeGb: {description: "The size of JVM heap for assembly and variant calling", category: "advanced"}
        nonJvmMemoryGb: {description: "The amount of memory in Gb to be requested besides JVM memory.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        vcf: {description: "VCF file including variant allele fractions."}
        vcfIndex: {description: "Index of output VCF."}
        assembly: {description: "The GRIDSS assembly BAM."}
        assemblyIndex: {description: "Index of output BAM file."}
    }
}

task GridssAnnotateVcfRepeatmasker {
    input {
        File gridssVcf
        File gridssVcfIndex
        String outputPath = "./gridss.repeatmasker_annotated.vcf.gz"

        String memory = "25G"
        Int threads = 8
        String dockerImage = "quay.io/biowdl/gridss:2.12.2"
        Int timeMinutes = 1440
    }

    command {
        gridss_annotate_vcf_repeatmasker \
        --output ~{outputPath} \
        --jar /usr/local/share/gridss-2.12.2-0/gridss.jar \
        -w . \
        -t ~{threads} \
        ~{gridssVcf}
    }

    output {
        File annotatedVcf = outputPath
        File annotatedVcfIndex = "~{outputPath}.tbi"
    }

    runtime {
        cpu: threads
        memory: memory
        time_minutes: timeMinutes # !UnknownRuntimeKey
        docker: dockerImage
    }

    parameter_meta {
        gridssVcf: {description: "The GRIDSS output.", category: "required"}
        gridssVcfIndex: {description: "The index for the GRIDSS output.", category: "required"}
        outputPath: {description: "The path the output should be written to.", category: "common"}
        threads: {description: "The number of the threads to use.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}

task Virusbreakend {
    input {
        File bam
        File bamIndex
        File referenceFasta
        File referenceFastaFai
        File referenceFastaDict
        File referenceImg
        File virusbreakendDB
        String outputPath = "./virusbreakend.vcf"

        String memory = "75G"
        Int threads = 12
        String dockerImage = "quay.io/biowdl/gridss:2.12.2"
        Int timeMinutes = 320
    }

    command {
        set -e
        mkdir virusbreakenddb
        tar -xzvf ~{virusbreakendDB} -C virusbreakenddb --strip-components 1
        virusbreakend \
        --output ~{outputPath} \
        --workingdir . \
        --reference ~{referenceFasta} \
        --db virusbreakenddb \
        --jar /usr/local/share/gridss-2.12.2-0/gridss.jar \
        -t ~{threads} \
        ~{bam}
    }

    output {
        File vcf = outputPath
        File summary = "~{outputPath}.summary.tsv"
    }

    runtime {
        cpu: threads
        memory: memory
        time_minutes: timeMinutes # !UnknownRuntimeKey
        docker: dockerImage
    }

    parameter_meta {
        bam: {description: "A BAM file.", category: "required"}
        bamIndex: {description: "The index for the BAM file.", category: "required"}
        referenceFasta: {description: "The fasta of the reference genome.", category: "required"}
        referenceImg: {description: "The BWA index image (generated with GATK BwaMemIndexImageCreator) of the reference.", category: "required"}
        virusbreakendDB: {description: "A .tar.gz containing the virusbreakend database.", category: "required"}
        outputPath: {description: "The path the output should be written to.", category: "common"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        threads: {description: "The number of the threads to use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}
