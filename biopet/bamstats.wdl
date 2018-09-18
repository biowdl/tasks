version 1.0

# Copyright Sequencing Analysis Support Core - Leiden University Medical Center 2018

import "../common.wdl" as common

task Generate {
    input {
        String? preCommand
        File? toolJar
        IndexedBamFile bam
        File? bedFile
        Boolean scatterMode = false
        Boolean onlyUnmapped = false
        Boolean tsvOutputs = false
        String outputDir
        Reference? reference
        Int memory = 4
        Float memoryMultiplier = 2.0
    }

    String toolCommand = if defined(toolJar)
        then "java -Xmx" + memory + "G -jar " + toolJar
        else "biopet-bamstats -Xmx" + memory + "G"

    String refArg = if (defined(reference)) then "--reference " + select_first([reference]).fasta else ""

    command {
        set -e -o pipefail
        ~{preCommand}
        mkdir -p ~{outputDir}
        ~{toolCommand} Generate \
        --bam ~{bam.file} \
        ~{"--bedFile " + bedFile} \
        ~{refArg} \
        ~{true="--onlyUnmapped" false="" onlyUnmapped} \
        ~{true="--scatterMode" false="" scatterMode} \
        ~{true="--tsvOutputs" false="" tsvOutputs} \
        --outputDir ~{outputDir}
    }

    output {
        File json = outputDir + "/bamstats.json"
        File summaryJson = outputDir + "/bamstats.summary.json"
    }

    runtime {
        memory: ceil(memory * memoryMultiplier)
    }
}