version 1.0

import "common.wdl" as common

task Germline {
    input {
        String runDir = "./strelka_run"
        Array[File]+ bams
        Array[File]+ indexes
        File referenceFasta
        File referenceFastaFai
        File? callRegions
        File? callRegionsIndex
        Boolean exome = false
        Boolean rna = false

        Int cores = 1
        Int memoryGb = 4
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
        docker: dockerImage
        cpu: cores
        memory: "~{memoryGb}G"
    }

    parameter_meta {
        runDir: {description: "The directory to use as run/output directory.", category: "common"}
        bams: {description: "The input BAM files.", category: "required"}
        indexes: {description: "The indexes for the input BAM files.", category: "required"}
        referenceFasta: {description: "The reference fasta file which was also used for mapping.", category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        callRegions: {description: "The bed file which indicates the regions to operate on.", category: "common"}
        callRegionsIndex: {description: "The index of the bed file which indicates the regions to operate on.", category: "common"}
        exome: {description: "Whether or not the data is from exome sequencing.", category: "common"}
        rna: {description: "Whether or not the data is from RNA sequencing.", category: "common"}

        cores: {description: "The number of cores to use.", category: "advanced"}
        memoryGb: {description: "The amount of memory this job will use in Gigabytes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
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
        File? callRegions
        File? callRegionsIndex
        File? indelCandidatesVcf
        File? indelCandidatesVcfIndex
        Boolean exome = false

        Int cores = 1
        Int memoryGb = 4
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
        docker: dockerImage
        cpu: cores
        memory: "~{memoryGb}G"
    }

    parameter_meta {
        runDir: {description: "The directory to use as run/output directory.", category: "common"}
        normalBam: {description: "The normal/control sample's BAM file.", category: "required"}
        normalBamIndex: {description: "The index for the normal/control sample's BAM file.", category: "required"}
        tumorBam: {description: "The tumor/case sample's BAM file.", category: "required"}
        tumorBamIndex: {description: "The index for the tumor/case sample's BAM file.", category: "required"}
        referenceFasta: {description: "The reference fasta file which was also used for mapping.", category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        callRegions: {description: "The bed file which indicates the regions to operate on.", category: "common"}
        callRegionsIndex: {description: "The index of the bed file which indicates the regions to operate on.", category: "common"}
        indelCandidatesVcf: {description: "An indel candidates VCF file from manta.", category: "advanced"}
        indelCandidatesVcfIndex: {description: "The index for the indel candidates VCF file.", category: "advanced"}
        exome: {description: "Whether or not the data is from exome sequencing.", category: "common"}

        cores: {description: "The number of cores to use.", category: "advanced"}
        memoryGb: {description: "The amount of memory this job will use in Gigabytes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }

    meta {
        WDL_AID: {
            exclude: ["doNotDefineThis"]
        }
    }
}