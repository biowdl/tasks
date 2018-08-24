version 1.0

task ConfigureGermline {
    input {
        String? preCommand
        String? installDir
        String runDir
        Array[File]+ bams
        Array[File]+ indexes
        File refFasta
        File refFastaIndex
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
        --ref ~{refFasta} \
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
        File normalBam
        File normalIndex
        File tumorBam
        File tumorIndex
        File refFasta
        File refFastaIndex
        File? callRegions
        File? callRegionsIndex
        File? indelCandidates
        File? indelCandidatesIndex
        Boolean exome = false
    }

    String toolCommand = if defined(installDir)
        then installDir + "bin/configureStrelkaSomaticWorkflow.py"
        else "configureStrelkaSomaticWorkflow.py"

    command {
        set -e -o pipefail
        ~{preCommand}
        ~{toolCommand} \
        --normalBam ~{normalBam} \
        --tumorBam ~{tumorBam} \
        --ref ~{refFasta} \
        --runDir ~{runDir} \
        ~{"--callRegions " + callRegions} \
        ~{"--indelCandidates " + indelCandidates} \
        ~{true="--exome" false="" exome} \
    }

    output {
        String runDirectory = runDir
    }
}

task Run {
    input {
        String? preCommand
        String runDir
        Int cores = 1
        Int memory = 4
        Boolean somatic = true
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