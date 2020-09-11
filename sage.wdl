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
        File tumorBamIndex
        String? normalName
        File? normalBam
        File? normalBamIndex
        File referenceFasta
        File referenceFastaDict
        File referenceFastaFai
        File hotspots
        File panelBed
        File highConfidenceBed
        Boolean hg38 = false
        String outputPath = "./sage.vcf.gz"

        Int threads = 2
        String javaXmx = "32G"
        String memory = "33G"
        Int timeMinutes = 1 + ceil(size(select_all([tumorBam, normalBam]), "G") * 10 / threads) #FIXME make sure this is enough
        String dockerImage = "quay.io/biocontainers/hmftools-sage:2.2--2"
    }

    command {
        java -Xmx~{javaXmx} -XX:ParallelGCThreads=1 \
        -cp /usr/local/share/hmftools-sage-2.2-2/sage.jar \
        com.hartwig.hmftools.sage.SageApplication \
        -tumor ~{tumorName} \
        -tumor_bam ~{tumorBam} \
        ~{"-reference " + normalName} \
        ~{"-reference_bam " + normalBam} \
        -ref_genome ~{referenceFasta} \
        -hotspots ~{hotspots} \
        -panel_bed ~{panelBed} \
        -high_confidence_bed ~{highConfidenceBed} \
        -assembly ~{true="hg38" false="hg19" hg38} \
        -threads ~{threads} \
        -out ~{outputPath}
    }

    output {
        File outputVcf = outputPath
        File outputVcfIndex = outputPath + ".tbi"
        # There is some plots as well, but in the current container the labels in the plots are just series of `□`s.
        # This seems to be a systemic issue with R generated plots in biocontainers...
    }

    runtime {
        time_minutes: timeMinutes # !UnknownRuntimeKey
        cpu: threads
        docker: dockerImage
        memory: memory
    }

    parameter_meta {
        tumorName: {description: "The name of the tumor sample.", category: "required"}
        tumorBam: {description: "The BAM file for the tumor sample.", category: "required"}
        tumorBamIndex: {description: "The index of the BAM file for the tumor sample.", category: "required"}
        normalName: {description: "The name of the normal/reference sample.", category: "common"}
        normalBam: {description: "The BAM file for the normal sample.", category: "common"}
        normalBamIndex: {description: "The index of the BAM file for the normal sample.", category: "common"}
        referenceFasta: {description: "The reference fasta file.", category: "required"}
        referenceFastaDict: {description: "The sequence dictionary associated with the reference fasta file.",
                             category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        hotspots: {description: "A vcf file with hotspot variant sites.", category: "required"}
        panelBed: {description: "A bed file describing coding regions to search for in frame indels.", category: "required"}
        highConfidenceBed: {description: "A bed files describing high confidence mapping regions.", category: "required"}

        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.",
                  category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}
