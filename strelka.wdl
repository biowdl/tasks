version 1.0

import "common.wdl" as common

task ConfigureGermline {
    input {
        String? preCommand
        String? installDir
        String runDir
        Array[File]+ bams
        Array[File]+ indexes
        Reference reference
        File? callRegions
        File? callRegionsIndex
        Boolean exome = false
        Boolean rna = false
    }

    String toolCommand = if defined(installDir)
        then installDir + "bin/configureStrelkaGermlineWorkflow.py"
        else "configureStrelkaGermlineWorkflow.py"

    command {
        set -e -o pipefail
        ~{preCommand}
        ~{toolCommand} \
        --bam ~{sep=" --bam " bams} \
        --ref ~{reference.fasta} \
        --runDir ~{runDir} \
        ~{"--callRegions " + callRegions} \
        ~{true="--exome" false="" exome} \
        ~{true="--rna" false="" rna}
    }

    output {
        String runDirectory = runDir
    }
}

task ConfigureSomatic {
    input {
        String? preCommand
        String? installDir
        String runDir
        IndexedBamFile normalBam
        IndexedBamFile tumorBam
        Reference reference
        File? callRegions
        File? callRegionsIndex
        IndexedVcfFile? indelCandidates
        Boolean exome = false
    }

    String toolCommand = if defined(installDir)
        then installDir + "bin/configureStrelkaSomaticWorkflow.py"
        else "configureStrelkaSomaticWorkflow.py"

    String indelCandidatesArg = if (defined(indelCandidates)) then "--indelCandidates " + select_first([indelCandidates]).file else ""

    command {
        set -e -o pipefail
        ~{preCommand}
        ~{toolCommand} \
        --normalBam ~{normalBam.file} \
        --tumorBam ~{tumorBam.file} \
        --ref ~{reference.fasta} \
        --runDir ~{runDir} \
        ~{"--callRegions " + callRegions} \
        ~{indelCandidatesArg} \
        ~{true="--exome" false="" exome} \
    }

    output {
        String runDirectory = runDir
    }
}

task Run {
    input {
        String runDir
        Int cores = 1
        Int memory = 4
        Boolean somatic = true
        #FIXME: This task does not have input files
    }

    command {
        ~{runDir}/runWorkflow.py \
        -m local \
        -j ~{cores} \
        -g ~{memory}
    }

    output {
        File? indelsVcf = runDir + "/results/variants/somatic.indels.vcf.gz"
        File? indelsIndex = runDir + "/results/variants/somatic.indels.vcf.gz.tbi"
        File variants = if somatic
            then runDir + "/results/variants/somatic.snvs.vcf.gz"
            else runDir + "/results/variants/variants.vcf.gz"
        File variantsIndex = if somatic
            then runDir + "/results/variants/somatic.snvs.vcf.gz.tbi"
            else runDir + "/results/variants/variants.vcf.gz.tbi"
    }

    runtime {
        cpu: cores
        memory: memory
    }
}