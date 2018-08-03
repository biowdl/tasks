version 1.0

task Germline {
    input {
        String? preCommand
        String? installDir
        String runDir
        Array[File]+ bams
        Array[File]+ indexes
        File refFasta
        File? callRegions
        File? callRegionsIndex
        Boolean exome = false
        Boolean rna = false

        Int cores = 1
        Int memory = 4
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
        ~{"--callRegions" + callRegions} \
        ~{true="--exome" false="" exome} \
        ~{true="--rna" false="" rna}

        ~{runDir}/runWorkflow.py \
        -m local \
        -J ~{cores} \
        -g ~{memory}
    }

    output {
        File variants = runDir + "/results/variants.vcf.gz"
        File variantsIndex = runDir + "/results/variants.vcf.gz.tbi"
    }

    runtime {
        cpu: cores
        memory: memory
    }
}


task Somatic {
    input {
        String? preCommand
        String? installDir
        String runDir
        File normalBam
        File normalIndex
        File tumorBam
        File tumorIndex
        File refFasta
        File? callRegions
        File? callRegionsIndex
        Boolean exome = false

        Int cores = 1
        Int memory = 4
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
        ~{"--callRegions" + callRegions} \
        ~{true="--exome" false="" exome} \

        ~{runDir}/runWorkflow.py \
        -m local \
        -J ~{cores} \
        -g ~{memory}
    }

    output {
        File indelsVcf = runDir + "/results/variants/somatic.indels.vcf.gz"
        File indelsIndex = runDir + "/results/variants/somatic.indels.vcf.gz.tbi"
        File snvVcf = runDir + "/results/variants/somatic.snvs.vcf.gz"
        File snvIndex = runDir + "/results/variants/somatic.snvs.vcf.gz.tbi"
    }

    runtime {
        cpu: cores
        memory: memory
    }
}