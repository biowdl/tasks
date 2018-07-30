version 1.0

task VarDict {
    input {
        String? installDir

        File tumorBam
        File normalBam
        File refFasta
        File bedFile
        String tumorSampleName
        String normalSampleName
        String outputVcf

        Int chromosomeColumn = 1
        Int startColumn = 2
        Int endColumn = 3
        Int geneColumn = 4

        String? preCommand
    }

    String toolCommand = if defined(installDir)
        then installDir + "/VarDict"
        else "vardict"

    command {
        set -e -o pipefail
        ~{preCommand}
        ~{toolCommand} \
        -G ~{refFasta} \
        -N ~{tumorSampleName} \
        -b "~{tumorBam}|~{normalBam}" \
        -c ~{chromosomeColumn} \
        -S ~{startColumn} \
        -E ~{endColumn} \
        -g ~{geneColumn} \
        ~{bedFile} | \
        ~{installDir + "/"}testsomatic.R | \
        ~{installDir + "/"}var2vcf_paired.pl \
        -N "~{tumorSampleName}|~{normalSampleName}" \
        > ~{outputVcf}
    }

    output {
        File vcfFile = outputVcf
    }
}