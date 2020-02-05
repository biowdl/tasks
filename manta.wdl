version 1.0

import "common.wdl"

task Somatic {
    input {
        File tumorBam
        File tumorBamIndex
        File? normalBam
        File? normalBamIndex
        File referenceFasta
        File referenceFastaFai
        String runDir = "./manta_run"
        File? callRegions
        File? callRegionsIndex
        Boolean exome = false

        Int cores = 1
        Int memory = 4
        String dockerImage = "quay.io/biocontainers/manta:1.4.0--py27_1"

    }

    command {
        configManta.py \
        ~{"--normalBam " + normalBam} \
        ~{"--tumorBam " + tumorBam} \
        --referenceFasta ~{referenceFasta} \
        ~{"--callRegions " + callRegions} \
        --runDir ~{runDir} \
        ~{true="--exome" false="" exome}

        ~{runDir}/runWorkflow.py \
        -m local \
        -j ~{cores} \
        -g ~{memory}
    }

    output {
        File candidateSmallIndelsVcf = runDir + "/results/variants/candidateSmallIndels.vcf.gz"
        File candidateSmallIndelsVcfIndex = runDir + "/results/variants/candidateSmallIndels.vcf.gz.tbi"
        File candidateSVVcf = runDir + "/results/variants/candidateSV.vcf.gz"
        File candidatSVVcfIndex = runDir + "/results/variants/candidateSV.vcf.gz.tbi"
        File tumorSVVcf = if defined(normalBam)
                          then runDir + "/results/variants/somaticSV.vcf.gz"
                          else runDir + "/results/variants/tumorSV.vcf.gz"
        File tumorSVVcfIndex = if defined(normalBam)
                               then runDir + "/results/variants/somaticSV.vcf.gz.tbi"
                               else runDir + "/results/variants/tumorSV.vcf.gz.tbi"
        File? diploidSV = runDir + "/results/variants/diploidSV.vcf.gz"
        File? diploidSVindex = runDir + "/results/variants/diploidSV.vcf.gz.tbi"
    }

    runtime {
        cpu: cores
        memory: memory
        docker: dockerImage
    }
}

task Germline {
    input {
        File bamFile
        File referenceFasta
        File referenceFastaFai
        String runDir
        File? callRegions
        File? callRegionsIndex
        Boolean exome = false
        
        Int cores = 1
        Int memory = 4
        String dockerTag = "1.4.0--py27_1"
    }

    command {
        set -e
        configManta.py \
        ~{"--normalBam " + bamFile} \
        --referenceFasta ~{referenceFasta} \
        ~{"--callRegions " + callRegions} \
        --runDir ~{runDir} \
        ~{true="--exome" false="" exome}
        
        ~{runDir}/runWorkflow.py \
        -m local \
        -j ~{cores} \
        -g ~{memory}
    }

    output {
        IndexedVcfFile candidateSmallIndels = object {
            file: runDir + "/results/variants/candidateSmallIndels.vcf.gz",
            index: runDir + "/results/variants/candidateSmallIndels.vcf.gz.tbi"
        }
        IndexedVcfFile candidateSV = object {
            file: runDir + "/results/variants/candidateSV.vcf.gz",
            index: runDir + "/results/variants/candidateSV.vcf.gz.tbi"
        }
        IndexedVcfFile diploidSV = object {
            file: runDir + "/results/variants/diploidSV.vcf.gz",
            index: runDir + "/results/variants/diploidSV.vcf.gz.tbi"
        }
    }
    
    runtime {
        cpu: cores
        memory: memory
        docker: "quay.io/biocontainers/manta:" + dockerTag
    }
}
    
