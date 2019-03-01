version 1.0

import "common.wdl" as common

task Germline {
    input {
        String runDir
        Array[File]+ bams
        Array[File]+ indexes
        Reference reference
        File? callRegions
        File? callRegionsIndex
        Boolean exome = false
        Boolean rna = false

        Int cores = 1
        Int memory = 4
        String dockerTag = "2.9.7--0"
    }

    command {
        set -e
        configureStrelkaGermlineWorkflow.py \
        --bam ~{sep=" --bam " bams} \
        --ref ~{reference.fasta} \
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
        docker: "quay.io/biocontainers/strelka:" + dockerTag
        cpu: cores
        memory: memory
    }
}

task Somatic {
    input {
        String runDir
        IndexedBamFile normalBam
        IndexedBamFile tumorBam
        Reference reference
        File? callRegions
        File? callRegionsIndex
        IndexedVcfFile? indelCandidates
        Boolean exome = false

        Int cores = 1
        Int memory = 4
        String dockerTag = "2.9.7--0"

        File? doNotDefineThis #FIXME
    }

    File? indelCandidatesFile = if (defined(indelCandidates))
        then select_first([indelCandidates]).file
        else doNotDefineThis

    command {
        set -e
        configureStrelkaSomaticWorkflow.py \
        --normalBam ~{normalBam.file} \
        --tumorBam ~{tumorBam.file} \
        --ref ~{reference.fasta} \
        --runDir ~{runDir} \
        ~{"--callRegions " + callRegions} \
        ~{"--indelCandidates " + indelCandidatesFile} \
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
        docker: "quay.io/biocontainers/strelka:" + dockerTag
        cpu: cores
        memory: memory
    }
}