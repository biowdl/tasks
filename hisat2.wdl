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

    String bamIndexPath = sub(outputBam, "\.bam$", ".bai")

    command {
        set -e -o pipefail
        mkdir -p "$(dirname ~{outputBam})"
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

    parameter_meta {
        indexFiles: {description: "The hisat2 index files.", category: "required"}
        inputR1: {description: "The first-/single-end FastQ file.", category: "required"}
        inputR2: {description: "The second-end FastQ file.", category: "common"}
        outputBam: {description: "The location the output BAM file should be written to.", category: "required"}
        sample: {description: "The sample id.", category: "required"}
        library: {description: "The library id.", category: "required"}
        readgroup: {description: "The readgroup id.", category: "required"}
        platform: {description: "The platform used for sequencing.", category: "advanced"}
        downstreamTranscriptomeAssembly: {description: "Equivalent to hisat2's `--dta` flag.", category: "advanced"}
        threads: {description: "The number of threads to use.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}