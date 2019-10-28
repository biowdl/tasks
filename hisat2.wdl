version 1.0

task Hisat2 {
    input {
        Array[File]+ indexFiles
        File inputR1
        File? inputR2
        String outputBam
        String sample
        String library
        String readgroup
        String platform = "illumina"
        Boolean downstreamTranscriptomeAssembly = true

        Int threads = 1
        String memory = "48G"
        # quay.io/biocontainers/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1
        # is a combination of hisat2 and samtools
        # hisat2=2.1.0, samtools=1.8
        String dockerImage = "quay.io/biocontainers/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:2388ff67fc407dad75774291ca5038f40cac4be0-0"
    }

    # Select_first is needed, otherwise womtool validate fails.
    String bamIndexPath = sub(outputBam, "\.bam$", ".bai")

    command {
        set -e -o pipefail
        mkdir -p $(dirname ~{outputBam})
        hisat2 \
        -p ~{threads} \
        -x ~{sub(indexFiles[0], "\.[0-9]\.ht2", "")} \
        ~{true="-1" false="-U" defined(inputR2)} ~{inputR1} \
        ~{"-2" + inputR2} \
        --rg-id ~{readgroup} \
        --rg 'SM:~{sample}' \
        --rg 'LB:~{library}' \
        --rg 'PL:~{platform}' \
        ~{true="--dta" false="" downstreamTranscriptomeAssembly} \
        | samtools sort > ~{outputBam}
        samtools index ~{outputBam} ~{bamIndexPath}
    }

    output {
        File bamFile = outputBam
        File bamIndex = bamIndexPath
    }

    runtime {
        memory: memory
        cpu: threads + 1
        docker: dockerImage
    }
}