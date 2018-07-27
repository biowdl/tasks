version 1.0

task Somatic {
    input {
        String? preCommand
        String? installDir
        String runDir
        File normalBam
        File tumorBam
        File refFasta

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
        --runDir ~{runDir}

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