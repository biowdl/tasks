version 1.0

# Copyright (c) 2025 Leiden University Medical Center
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

task Pileup {
    input {
        String dockerImage = "quay.io/biocontainers/ont-modkit:0.4.2--hcdda2d0_0"
        File bam
        File bamIndex
        String outputBed = "output.bed"
        File referenceFasta
        File referenceFastaFai

        Int? intervalSize
        File? includeBed

        Boolean cpg = false
        Boolean combineMods = false
        Boolean combineStrands = false
        Boolean bedgraph = false
        String? ignore
        String logFilePath = "modkit.log"

        Int threads = 8
        String memory = "48GiB"
        Int timeMinutes = 4320  # 3 Days

    }

    command <<<
        set -e
        mkdir -p $(dirname ~{outputBed})
        mkdir -p $(dirname ~{logFilePath})
        modkit pileup \
        --threads ~{threads} \
        ~{"--interval-size " + intervalSize} \
        ~{"--include-bed " + includeBed} \
        ~{"--ignore " + ignore} \
        --ref ~{referenceFasta} \
        ~{true="--cpg" false="" cpg} \
        ~{true="--combine-mods" false="" combineMods} \
        ~{true="--combine-strands" false="" combineStrands} \
        ~{true="--bedgraph" false="" bedgraph} \
        --log-filepath ~{logFilePath} \
        ~{bam} \
        ~{outputBed} 
    >>>

    output {
        File out = outputBed
        File logFile = logFilePath
    }

    runtime {
        docker: dockerImage
        cpu: threads
        memory: memory
        time_minutes: timeMinutes
    }
}