version 1.0

import "common.wdl" as common

task Germline {
    input {
        String runDir = "."
        Array[File]+ bams
        Array[File]+ indexes
        File referenceFasta
        File? callRegions
        File? callRegionsIndex
        Boolean exome = false
        Boolean rna = false

        Int cores = 1
        Int memory = 4
        String dockerImage = "quay.io/biocontainers/strelka:2.9.7--0"
    }

    command {
        set -e
        mkdir -p ~{runDir}
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
        -g ~{memory}
    }

    output {
        File variants = runDir + "/results/variants/variants.vcf.gz"
        File variantsIndex = runDir + "/results/variants/variants.vcf.gz.tbi"
    }

    runtime {
        docker: dockerImage
        cpu: cores
        memory: memory
    }
}

task Somatic {
    input {
        String runDir = "."
        File normalBam
        File normalBamIndex
        File tumorBam
        File tumorBamIndex
        File referenceFasta
        File? callRegions
        File? callRegionsIndex
        File? indelCandidatesVcf
        File? indelCandidatesVcfIndex
        Boolean exome = false

        Int cores = 1
        Int memory = 4
        String dockerImage = "quay.io/biocontainers/strelka:2.9.7--0"

        File? doNotDefineThis #FIXME
    }

    command {
        set -e
        mkdir -p ~{runDir}
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
        -g ~{memory}
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
        memory: memory
    }
}