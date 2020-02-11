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

        String memory = "16G"
        String javaXmx = "8G"
    }

    File referenceFasta = if defined(reference) then select_first([reference]).fasta else ""

    String toolCommand = if defined(toolJar)
        then "java -Xmx~{javaXmx} -jar " + toolJar
        else "biopet-bamstats -Xmx~{javaXmx}"

    command {
        set -e -o pipefail
        ~{preCommand}
        mkdir -p ~{outputDir}
        ~{toolCommand} Generate \
        --bam ~{bam.file} \
        ~{"--bedFile " + bedFile} \
        ~{true="--reference" false="" defined(reference)} ~{referenceFasta} \
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
        memory: memory
    }
}