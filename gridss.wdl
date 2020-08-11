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

        Int threads = 1
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
        memory: "32G"
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
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}
    }
}