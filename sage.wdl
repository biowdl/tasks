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

task SageHotspot {
    input {
        String tumorName
        File tumorBam
        File tumorBai
        String? normalName
        File? normalBam
        File? normalBai
        File referenceFasta
        File referenceFastaDict
        File referenceFastaFai
        File knownHotspots
        File codingRegsions
        String outputPath = "./sage_hotspot.vcf.gz"

        Int timeMinutes = 60 #FIXME I've no idea how long this takes...
        String javaXmx = "32G"
        String memory = "33G"
        String dockerImage = "quay.io/biocontainers/hmftools-sage:2.2--0"
    }

    command {
        SAGE -Xmx~{javaXmx} \
        -tumor ~{tumorName} \
        -tumor_bam ~{tumorBam} \
        ~{"-reference " + normalName} \
        ~{"-reference_bam " + normalBam} \
        -ref_genome ~{referenceFasta} \
        -known_hotspots ~{knownHotspots} \
        -coding_regions ~{codingRegsions} \
        -out ~{outputPath}
    }

    output {
        File outputVcf = outputPath
    }

    runtime {
        time_minutes: timeMinutes
        docker: dockerImage
        memory: memory
    }

    parameter_meta {
        tumorName: {description: "The name of the tumor sample.", category: "required"}
        tumorBam: {description: "The BAM file for the tumor sample.", category: "required"}
        tumorBai: {description: "The index of the BAM file for the tumor sample.", category: "required"}
        normalName: {description: "The name of the normal/reference sample.", category: "common"}
        normalBam: {description: "The BAM file for the normal sample.", category: "common"}
        normalBai: {description: "The index of the BAM file for the normal sample.", category: "common"}
        referenceFasta: {description: "The reference fasta file.", category: "required"}
        referenceFastaDict: {description: "The sequence dictionary associated with the reference fasta file.",
                             category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        knownHotspots: {description: "A TSV file with hotspot variant sites.", category: "required"}
        codingRegsions: {description: "A bed file describing coding regions to search for inframe indels.", category: "required"}

        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.",
                  category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}

task Sage {
    input {
        String tumorName
        File tumorBam
        String? normalName
        File? normalBam
        File referenceFasta
        File referenceFastaFai
        File referenceFastaDict
        File hotspots
        File panelBed
        File highConfidenceBed
        String assembly = "hg38"
        String outputPath = "./sage.vcf.gz"

        Int timeMinutes = 60 #FIXME I've no idea how long this takes...
        String javaXmx = "32G"
        String memory = "33G"
        Int threads = 2
        String dockerImage = "quay.io/biocontainers/hmftools-sage:2.2--0"
    }

    command {
        java -Xmx~{javaXmx} \
        -cp /usr/local/share/hmftools-sage-2.2-0/sage.jar \
        com.hartwig.hmftools.sage.SageApplication \
        -tumor ~{tumorName} \
        -tumor_bam ~{tumorBam} \
        ~{"-reference " + normalName} \
        ~{"-reference_bam " + normalBam} \
        -ref_genome ~{referenceFasta} \
        -hotspots ~{hotspots} \
        -panel_bed ~{panelBed} \
        -high_confidence_bed ~{highConfidenceBed} \
        -assembly ~{assembly} \
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
        normalBai: {description: "The index of the BAM file for the normal sample.", category: "common"}
        referenceFasta: {description: "The reference fasta file.", category: "required"}
        referenceFastaDict: {description: "The sequence dictionary associated with the reference fasta file.",
                             category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        hotspots: {description: "A VCF file containg hotspot variant sites.", category: "required"}
        panelBed: {description: "A bed file containing a panel of genes of intrest.", category: "required"}
        highConfidenceBed: {description: "A bed file containing high confidence regions.", category: "required"}
        assembly: {description: "The genome assembly used, either \"hg19\" or \"hg38\".", category: "common"}
        outputPath: {description: "The path to write the output VCF to.", category: "common"}

        threads: {description: "The number of threads to be used.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.",
                  category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}