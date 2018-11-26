version 1.0

import "common.wdl"

task ConfigureSomatic {
    input {
        IndexedBamFile tumorBam
        IndexedBamFile? normalBam
        Reference reference
        String runDir
        File? callRegions
        File? callRegionsIndex
        Boolean exome = false
        String? preCommand
        String? installDir
    }

    String toolCommand = if defined(installDir)
        then installDir + "bin/configMata.py"
        else "configManta.py"

    String normalArg = if (defined(normalBam)) then "--normalBam " + select_first([normalBam]).file else ""

    command {
        set -e -o pipefail
        ~{preCommand}
        ~{toolCommand} \
        ~{normalArg} \
        ~{"--tumorBam " + tumorBam.file} \
        --referenceFasta ~{reference.fasta} \
        ~{"--callRegions " + callRegions} \
        --runDir ~{runDir} \
        ~{true="--exome" false="" exome}
    }

    output {
        String runDirectory = runDir
    }
}

task RunSomatic {
    input {
        String? preCommand
        String runDir
        Int cores = 1
        Int memory = 4
        Boolean paired = true
    }

    command {
        set -e -o pipefail
        ~{preCommand}
        ~{runDir}/runWorkflow.py \
        -m local \
        -j ~{cores} \
        -g ~{memory}
    }

    output {
        IndexedVcfFile condidateSmallIndels = object {
          file: runDir + "/results/variants/candidateSmallIndels.vcf.gz",
          index: runDir + "/results/variants/candidateSmallIndels.vcf.gz.tbi"
        }
        IndexedVcfFile candidateSV = object {
          file: runDir + "/results/variants/candidateSV.vcf.gz",
          index: runDir + "/results/variants/candidateSV.vcf.gz.tbi"
        }
        IndexedVcfFile tumorSV = if (paired)
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
    }
}