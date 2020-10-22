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

task GRIDSS {
    input {
        File tumorBam
        File tumorBai
        String tumorLabel
        File? normalBam
        File? normalBai
        String? normalLabel
        BwaIndex reference
        String outputPrefix = "gridss"

        Int jvmHeapSizeGb = 30
        Int threads = 2
        Int timeMinutes = ceil(1440 / threads) + 10
        String dockerImage = "quay.io/biocontainers/gridss:2.9.4--0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPrefix})"
        gridss \
        --reference ~{reference.fastaFile} \
        --output ~{outputPrefix}.vcf.gz \
        --assembly ~{outputPrefix}_assembly.bam \
        ~{"-t " + threads} \
        ~{"--jvmheap " + jvmHeapSizeGb + "G"} \
        --label ~{normalLabel}~{true="," false="" defined(normalLabel)}~{tumorLabel} \
        ~{normalBam} \
        ~{tumorBam}
        tabix -p vcf ~{outputPrefix}.vcf.gz
        samtools index ~{outputPrefix}_assembly.bam ~{outputPrefix}_assembly.bai
    }

    output {
        File vcf = outputPrefix + ".vcf.gz"
        File vcfIndex = outputPrefix + ".vcf.gz.tbi"
        File assembly = outputPrefix + "_assembly.bam"
        File assemblyIndex = outputPrefix + "_assembly.bai"
    }

    runtime {
        cpu: threads
        memory: "~{jvmHeapSizeGb + 1}G"
        time_minutes: timeMinutes # !UnknownRuntimeKey
        docker: dockerImage
    }

    parameter_meta {
        tumorBam: {description: "The input BAM file. This should be the tumor/case sample in case of a paired analysis.", category: "required"}
        tumorBai: {description: "The index for tumorBam.", category: "required"}
        tumorLabel: {description: "The name of the (tumor) sample.", category: "required"}
        normalBam: {description: "The BAM file for the normal/control sample.", category: "advanced"}
        normalBai: {description: "The index for normalBam.", category: "advanced"}
        normalLabel: {description: "The name of the normal sample.", category: "advanced"}
        reference: {description: "A BWA index, this should also include the fasta index file (.fai).", category: "required"}
        outputPrefix: {description: "The prefix for the output files. This may include parent directories.", category: "common"}

        threads: {description: "The number of the threads to use.", category: "advanced"}
        jvmHeapSizeGb: {description: "The size of JVM heap for assembly and variant calling",category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}
    }
}

task AnnotateInsertedSequence {
    input {
        File inputVcf
        String outputPath = "gridss.annotated.vcf.gz"
        File viralReference

        Int threads = 8
        String javaXmx = "8G"
        String memory = "9G"
        String dockerImage = "quay.io/biocontainers/gridss:2.9.4--0"
        Int timeMinutes = 1 + ceil(size(inputVcf, "G") * 2 / threads)
    }

    command {
        java -Xmx~{javaXmx} \
        -Dsamjdk.create_index=true \
        -Dsamjdk.use_async_io_read_samtools=true \
        -Dsamjdk.use_async_io_write_samtools=true \
        -Dsamjdk.use_async_io_write_tribble=true \
        -Dsamjdk.buffer_size=4194304 \
        -cp /usr/local/share/gridss-2.9.4-0/gridss.jar \
        gridss.AnnotateInsertedSequence \
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

        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.",
                  category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}