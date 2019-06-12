version 1.0

task Hisat2 {
    input {
        String indexBasename
        File inputR1
        File? inputR2
        String outputBam
        String sample
        String library
        String readgroup
        String platform = "illumina"
        Boolean downstreamTranscriptomeAssembly = true

        Int threads = 1
        Int memory = 48
        # hisat2=2.1.0, samtools=1.8
        String dockerTag = "2388ff67fc407dad75774291ca5038f40cac4be0-0"
    }

    Array[File] indexFiles = glob(indexBasename + "*")
    String indexPath = sub(indexFiles[0], "\.[0-9]\.ht2", "")

    command {
        set -e -o pipefail
        mkdir -p $(dirname ~{outputBam})
        hisat2 \
        -p ~{threads} \
        -x ~{indexPath} \
        ~{true="-1" false="-U" defined(inputR2)} ~{inputR1} \
        ~{"-2" + inputR2} \
        --rg-id ~{readgroup} \
        --rg 'SM:~{sample}' \
        --rg 'LB:~{library}' \
        --rg 'PL:~{platform}' \
        ~{true="--dta" false="" downstreamTranscriptomeAssembly} \
        | samtools sort > ~{outputBam}
    }

    output {
        File bamFile = outputBam
    }

    runtime {
        memory: (memory / threads) + 1
        cpu: threads + 1
        # quay.io/biocontainers/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1
        # is a combination of hisat2 and samtools
        docker: "quay.io/biocontainers/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:" + dockerTag
    }
}