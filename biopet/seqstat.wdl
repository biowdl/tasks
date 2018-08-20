version 1.0

# Copyright Sequencing Analysis Support Core - Leiden University Medical Center 2018

task Generate {
    input {
        String? preCommand
        File? toolJar
        File fastqR1
        File? fastqR2
        String outputFile
        String sample = "sample"
        String library = "library"
        String readgroup = "readgroup"

        Int memory = 4
        Float memoryMultiplier = 2.0
    }

    String toolCommand = if defined(toolJar)
        then "java -Xmx" + memory + "G -jar " + toolJar
        else "biopet-seqstat -Xmx" + memory + "G"

    command {
        set -e -o pipefail
        ~{preCommand}
        mkdir -p $(dirname ~{outputFile})
        ~{toolCommand} Generate \
        --fastqR1 ~{fastqR1} \
        ~{"--fastqR2 " + fastqR2} \
        --output ~{outputFile} \
        ~{"--sample " + sample} \
        ~{"--library " + library } \
        ~{"--readgroup " + readgroup }
    }

    output {
        File json = outputFile
    }

    runtime {
        memory: ceil(memory * memoryMultiplier)
    }
}