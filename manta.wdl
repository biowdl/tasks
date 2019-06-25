version 1.0

import "common.wdl"

task Somatic {
    input {
        IndexedBamFile tumorBam
        IndexedBamFile? normalBam
        Reference reference
        String runDir
        File? callRegions
        File? callRegionsIndex
        Boolean exome = false

        Int cores = 1
        Int memory = 4
        String dockerImage = "quay.io/biocontainers/manta:1.4.0--py27_1"

        File? doNotDefineThis #FIXME
    }

    File? normalBamFile = if defined(normalBam)
        then select_first([normalBam]).file
        else doNotDefineThis

    command {
        set -e
        configManta.py \
        ~{"--normalBam " + normalBamFile} \
        ~{"--tumorBam " + tumorBam.file} \
        --referenceFasta ~{reference.fasta} \
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
        IndexedVcfFile tumorSV = if defined(normalBam)
            then object {
              file: runDir + "/results/variants/somaticSV.vcf.gz",
              index: runDir + "/results/variants/somaticSV.vcf.gz.tbi"
            }
            else object {
              file: runDir + "/results/variants/tumorSV.vcf.gz",
              index: runDir + "/results/variants/tumorSV.vcf.gz.tbi"
            }

        #FIXME: workaround for https://github.com/broadinstitute/cromwell/issues/4111
        File? diploidSV = runDir + "/results/variants/diploidSV.vcf.gz"
        File? diploidSVindex = runDir + "/results/variants/diploidSV.vcf.gz.tbi"
    }

    runtime {
        cpu: cores
        memory: memory
        docker: dockerImage
    }
}
