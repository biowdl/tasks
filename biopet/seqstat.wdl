version 1.0

# Copyright Sequencing Analysis Support Core - Leiden University Medical Center 2018

import "../common.wdl" as common

task Generate {
    input {
        String? preCommand
        File? toolJar
        FastqPair fastq
        String outputFile
        String sample
        String library
        String readgroup

        String memory = "10G"
        String javaXmx = "4G"
    }

    String toolCommand = if defined(toolJar)
        then "java -Xmx~{javaXmx} -jar " + toolJar
        else "biopet-seqstat -Xmx~{javaXmx}"

    command {
        set -e -o pipefail
        ~{preCommand}
        mkdir -p $(dirname ~{outputFile})
        ~{toolCommand} Generate \
        --fastqR1 ~{fastq.R1} \
        ~{"--fastqR2 " + fastq.R2} \
        --output ~{outputFile} \
        ~{"--sample " + sample} \
        ~{"--library " + library } \
        ~{"--readgroup " + readgroup }
    }

    output {
        File json = outputFile
    }

    runtime {
        memory: memory
    }
}