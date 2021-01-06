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

import "common.wdl" as common

task Germline {
    input {
        String runDir = "./strelka_run"
        Array[File]+ bams
        Array[File]+ indexes
        File referenceFasta
        File referenceFastaFai
        Boolean exome = false
        Boolean rna = false

        File? callRegions
        File? callRegionsIndex

        Int cores = 1
        Int memoryGb = 4
        Int timeMinutes = 90
        String dockerImage = "quay.io/biocontainers/strelka:2.9.7--0"
    }

    command {
        configureStrelkaGermlineWorkflow.py \
        --bam ~{sep=" --bam " bams} \
        --ref ~{referenceFasta} \
        --runDir ~{runDir} \
        ~{"--callRegions " + callRegions} \
        ~{true="--exome" false="" exome} \
        ~{true="--rna" false="" rna}

        ~{runDir}/runWorkflow.py \
        -m local \
        -j ~{cores} \
        -g ~{memoryGb}
    }

    output {
        File variants = runDir + "/results/variants/variants.vcf.gz"
        File variantsIndex = runDir + "/results/variants/variants.vcf.gz.tbi"
    }

    runtime {
        cpu: cores
        memory: "~{memoryGb}G"
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        runDir: {description: "The directory to use as run/output directory.", category: "common"}
        bams: {description: "The input BAM files.", category: "required"}
        indexes: {description: "The indexes for the input BAM files.", category: "required"}
        referenceFasta: {description: "The reference fasta file which was also used for mapping.", category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        exome: {description: "Whether or not the data is from exome sequencing.", category: "common"}
        rna: {description: "Whether or not the data is from RNA sequencing.", category: "common"}
        callRegions: {description: "The bed file which indicates the regions to operate on.", category: "common"}
        callRegionsIndex: {description: "The index of the bed file which indicates the regions to operate on.", category: "common"}
        cores: {description: "The number of cores to use.", category: "advanced"}
        memoryGb: {description: "The amount of memory this job will use in Gigabytes.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        variants: {description: "Output VCF file."}
        variantsIndex: {description: "Index of output VCF file."}
    }
}

task Somatic {
    input {
        String runDir = "./strelka_run"
        File normalBam
        File normalBamIndex
        File tumorBam
        File tumorBamIndex
        File referenceFasta
        File referenceFastaFai
        Boolean exome = false

        File? callRegions
        File? callRegionsIndex
        File? indelCandidatesVcf
        File? indelCandidatesVcfIndex

        Int cores = 1
        Int memoryGb = 4
        Int timeMinutes = 90
        String dockerImage = "quay.io/biocontainers/strelka:2.9.7--0"

        File? doNotDefineThis #FIXME
    }

    command {
        configureStrelkaSomaticWorkflow.py \
        --normalBam ~{normalBam} \
        --tumorBam ~{tumorBam} \
        --ref ~{referenceFasta} \
        --runDir ~{runDir} \
        ~{"--callRegions " + callRegions} \
        ~{"--indelCandidates " + indelCandidatesVcf} \
        ~{true="--exome" false="" exome}

        ~{runDir}/runWorkflow.py \
        -m local \
        -j ~{cores} \
        -g ~{memoryGb}
    }

    output {
        File indelsVcf = runDir + "/results/variants/somatic.indels.vcf.gz"
        File indelsIndex = runDir + "/results/variants/somatic.indels.vcf.gz.tbi"
        File variants = runDir + "/results/variants/somatic.snvs.vcf.gz"
        File variantsIndex = runDir + "/results/variants/somatic.snvs.vcf.gz.tbi"
    }

    runtime {
        cpu: cores
        memory: "~{memoryGb}G"
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        runDir: {description: "The directory to use as run/output directory.", category: "common"}
        normalBam: {description: "The normal/control sample's BAM file.", category: "required"}
        normalBamIndex: {description: "The index for the normal/control sample's BAM file.", category: "required"}
        tumorBam: {description: "The tumor/case sample's BAM file.", category: "required"}
        tumorBamIndex: {description: "The index for the tumor/case sample's BAM file.", category: "required"}
        referenceFasta: {description: "The reference fasta file which was also used for mapping.", category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        exome: {description: "Whether or not the data is from exome sequencing.", category: "common"}
        callRegions: {description: "The bed file which indicates the regions to operate on.", category: "common"}
        callRegionsIndex: {description: "The index of the bed file which indicates the regions to operate on.", category: "common"}
        indelCandidatesVcf: {description: "An indel candidates VCF file from manta.", category: "advanced"}
        indelCandidatesVcfIndex: {description: "The index for the indel candidates VCF file.", category: "advanced"}
        cores: {description: "The number of cores to use.", category: "advanced"}
        memoryGb: {description: "The amount of memory this job will use in Gigabytes.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        indelsVcf: {description: "VCF containing INDELS."}
        indelsIndex: {description: "Index of output `indelsVcf`."}
        variants: {description: "VCF containing variants."}
        variantsIndex: {description: "Index of output `variants`."}
    }

    meta {
        WDL_AID: {
            exclude: ["doNotDefineThis"]
        }
    }
}
