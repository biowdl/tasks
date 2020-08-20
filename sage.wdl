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

task Sage {
    input {
        String tumorName
        File tumorBam
        File tumorBai
        String? normalName
        File? normalBam
        File? normalBai
        String assembly
        File referenceFasta
        File hotspotVcf
        File panelBed
        File highConfidenceBed

        Int timeMinutes = 60 #FIXME I've no idea how long this takes...
        Int threads = 2
        String javaXmx = "32G"
        String dockerImage = "quay.io/biocontainers/hmftools-sage:2.2--0"
    }

    command {
        SAGE \
        -Xmx~{javaXmx} \
        -tumor ~{tumorName} \
        -tumor_bam ~{tumorBam} \
        ~{"-reference " + normalName} \
        ~{"-reference_bam " + normalBam} \
        -assembly ~{assembly} \
        -ref_genome ~{referenceFasta} \
        -hotspots ~{hotspotVcf} \
        -panel_bed ~{panelBed} \
        -high_confidence_bed ~{highConfidenceBed} \
        -threads ~{threads} \

        -out ~{outputPath}
    }

    output {
        File outputVcf = outputPath
    }

    runtime {
        time_minutes: timeMinutes
        cpu: threads
        docker: dockerImage
        memory: memory
    }

    parameter_meta {
        tumorName: {description: "The name of the tumor sample.", category: "required"}
        tumorBam: {description: "The BAM file for the tumor sample.", category: "required"}
        tumorBai: {description: "The index of the BAM file for the tumor sample.", category: "required"}
        normalName: {description: "The name of the normal/reference sample.", category: "common"}
        normalBam: {description: "The BAM file for the normal sample.", category: "common"}
        normalBam: {description: "The BAM file for the normal sample.", category: "common"}
        assembly: {description: "The assembly of the reference genomes, either hg19 or hg38.", category: "required"}
        referenceFasta: {description: "The reference fasta file.", category: "required"}
        hotspotVcf: {description: "A VCF file with hotspot variant sites.", category: "required"}
        panelBed: {description: "A bed file describing a panel of cancer related genes.", category: "required"}
        highConfidenceBed: {description: "A bed file describing high confidence regions.", category: "required"}

        threads: {description: "The number of threads to be used.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.",
                  category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}