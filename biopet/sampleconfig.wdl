version 1.0

import "../common.wdl" as common

task SampleConfig {
    input {
        File? toolJar
        String? preCommand
        Array[File]+ inputFiles
        String keyFilePath
        String? sample
        String? library
        String? readgroup
        String? jsonOutputPath
        String? tsvOutputPath

        String memory = "8G"
        String javaXmx = "16G"
    }

    String toolCommand = if defined(toolJar)
        then "java -Xmx~{javaXmx} -jar " + toolJar
        else "biopet-sampleconfig -Xmx~{javaXmx}"

    command {
        set -e -o pipefail
        ~{preCommand}
        mkdir -p . ~{"$(dirname " + jsonOutputPath + ")"} ~{"$(dirname " + tsvOutputPath + ")"}
        ~{toolCommand} \
        -i ~{sep="-i " inputFiles} \
        ~{"--sample " + sample} \
        ~{"--library " + library} \
        ~{"--readgroup " + readgroup} \
        ~{"--jsonOutput " + jsonOutputPath} \
        ~{"--tsvOutput " + tsvOutputPath} \
        > ~{keyFilePath}
    }

    output {
        File keysFile = keyFilePath
        File? jsonOutput = jsonOutputPath
        File? tsvOutput = tsvOutputPath
    }

    runtime {
        memory: memory
    }
}

task SampleConfigCromwellArrays {
    input {
        File? toolJar
        String? preCommand
        Array[File]+ inputFiles
        String outputPath

        String memory = "8G"
        String javaXmx = "4G"
    }

    String toolCommand = if defined(toolJar)
        then "java -Xmx~{javaXmx} -jar " + toolJar
        else "biopet-sampleconfig -Xmx~{javaXmx}"

    command {
        set -e -o pipefail
        ~{preCommand}
        mkdir -p $(dirname ~{outputPath})
        ~{toolCommand} CromwellArrays \
        -i ~{sep="-i " inputFiles} \
        ~{"-o " + outputPath}
    }

    output {
        File outputFile = outputPath
    }

    runtime {
        memory: memory
    }
}

task CaseControl {
    input {
        File? toolJar
        String? preCommand
        Array[File]+ inputFiles
        Array[File]+ inputIndexFiles
        Array[File]+ sampleConfigs
        String outputPath
        String controlTag = "control"

        String memory = "8G"
        String javaXmx = "4G"
    }

    String toolCommand = if defined(toolJar)
        then "java -Xmx~{javaXmx} -jar " + toolJar
        else "biopet-sampleconfig -Xmx~{javaXmx}"

    command {
        set -e -o pipefail
        ~{preCommand}
        mkdir -p $(dirname ~{outputPath})
        ~{toolCommand} CaseControl \
        -i ~{sep=" -i " inputFiles} \
        -s ~{sep=" -s " sampleConfigs} \
        ~{"-o " + outputPath} \
        ~{"--controlTag " + controlTag}
    }

    output {
        File outputFile = outputPath
        CaseControls caseControls = read_json(outputFile)
    }

    runtime {
        memory: memory
    }
}
