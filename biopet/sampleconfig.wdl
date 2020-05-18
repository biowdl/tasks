version 1.0

# Copyright (c) 2017 Leiden University Medical Center
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

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

        String memory = "17G"
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

        String memory = "5G"
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

        String memory = "5G"
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
