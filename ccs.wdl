version 1.0

# Copyright (c) 2020 Sequencing Analysis Support Core - Leiden University Medical Center
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

task CCS {
    input {
        Int minPasses = 3
        Int minLength = 10
        Int maxLength = 50000
        Boolean byStrand = false
        Float minReadQuality = 0.99
        String logLevel = "WARN"
        File subreadsFile
        String outputPrefix
        
        Int cores = 4
        String memory = "10G"
        String dockerImage = "quay.io/biocontainers/pbccs:4.2.0--0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPrefix})"
        ccs \
        ~{"--min-passes " + minPasses} \
        ~{"--min-length " + minLength} \
        ~{"--max-length " + maxLength} \
        ~{true="--by-strand" false="" byStrand} \
        ~{"--min-rq " + minReadQuality} \
        ~{"--log-level " + logLevel} \
        ~{"--num-threads " + cores} \
        ~{"--report-file " + outputPrefix + ".ccs.report.txt"} \
        ~{"--log-file " + outputPrefix + ".ccs.stderr.log"} \
        ~{subreadsFile}
        ~{outputPrefix + ".ccs.bam"}
    }

    output {
        File outputCCSfile = outputPrefix + ".ccs.bam"
        File outputCCSindexFile = outputPrefix + ".ccs.bam.pbi"
        File outputReportFile = outputPrefix + ".ccs.report.txt"
        File outputSTDERRfile = outputPrefix + ".ccs.stderr.log"
    }

    runtime {
        cpu: cores
        memory: memory
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        minPasses: {description: "Minimum number of full-length subreads required to generate CCS for a ZMW.", category: "advanced"}
        minLength: {description: "Minimum draft length before polishing.", category: "advanced"}
        maxLength: {description: "Maximum draft length before polishing.", category: "advanced"}
        byStrand: {description: "Generate a consensus for each strand.", category: "advanced"}
        minReadQuality: {description: "Minimum predicted accuracy in [0, 1].", category: "common"}
        logLevel: {description: "Set log level. Valid choices: (TRACE, DEBUG, INFO, WARN, FATAL).", category: "advanced"}
        subreadsFile: {description: "Subreads input file.", category: "required"}
        outputPrefix: {description: "Output directory path + output file prefix.", category: "required"}
        cores: {description: "The number of cores to be used.", category: "advanced"}
        memory: {description: "The amount of memory available to the job.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        outputCCSfile: {description: "Consensus reads output file."}
        outputCCSindexFile: {description: "Index of consensus reads output file."}
        outputReportFile: {description: "CCS results report file."}
        outputSTDERRfile: {description: "CCS STDERR log file."}
    }
}
