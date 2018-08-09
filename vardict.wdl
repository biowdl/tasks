version 1.0

task VarDict {
    input {
        String? installDir
        Boolean useJavaVersion = true

        String tumorSampleName
        File tumorBam
        File tumorIndex
        String? normalSampleName
        File? normalBam
        File? normalIndex
        File refFasta
        File refFastaIndex
        File bedFile
        String outputVcf

        Int chromosomeColumn = 1
        Int startColumn = 2
        Int endColumn = 3
        Int geneColumn = 4

        String? preCommand
        Int memory = 4
        Float memoryMultiplier = 2.0
    }

    String toolCommand = if defined(installDir)
        then installDir + "/VarDict"
        else if useJavaVersion
            then "vardict-java -Xmx${memory}"
            else "vardict"

    command {
        set -e -o pipefail
        ~{preCommand}
        ~{toolCommand} \
        -G ~{refFasta} \
        -N ~{tumorSampleName} \
        -b "~{tumorBam}~{"|" + normalBam}" \
        ~{true="-z" false="" defined(normalBam)} \
        -c ~{chromosomeColumn} \
        -S ~{startColumn} \
        -E ~{endColumn} \
        -g ~{geneColumn} \
        ~{bedFile} | \
        ~{installDir + "/"}~{true="testsomatic.R" false="teststrandbias.R" defined(normalBam)} | \
        ~{installDir + "/"}~{true="var2vcf_paired.pl"
            false="var2vcf_valid.pl" defined(normalBam)} \
        -N "~{tumorSampleName}~{"|" + normalSampleName}" \
        ~{true="" false="-E" defined(normalBam)} | \
        bgzip -c > ~{outputVcf}
    }

    output {
        File vcfFile = outputVcf
    }

    runtime {
        memory: ceil(memory * memoryMultiplier)
    }
}
