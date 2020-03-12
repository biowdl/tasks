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

task Refine {
    input {
        Int minPolyAlength = 20
        Boolean requirePolyA = false
        String logLevel = "WARN"
        File inputBamFile
        File primerFile
        String outputPrefix

        Int cores = 4
        String memory = "10G"
        String dockerImage = "quay.io/biocontainers/isoseq3:3.3.0--0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPrefix})"
        isoseq3 refine \
        ~{"--min-poly-length" + minPolyAlength} \
        ~{true="--require-polya" false="" requirePolyA} \
        ~{"--log-level " + logLevel} \
        ~{"--num-threads " + cores} \
        ~{"--log-file " + outputPrefix + ".flnc.stderr.log"} \
        ~{inputBamFile} \
        ~{primerFile} \
        ~{outputPrefix + ".flnc.bam"}
    }

    output {
        File outputFLfile = outputPrefix + ".flnc.bam"
        File outputFLindexFile = outputPrefix + ".flnc.bam.pbi"
        File outputSTDERRfile = outputPrefix + ".flnc.stderr.log"
        File outputConsensusReadsetFile = outputPrefix + ".consensusreadset.xml"
        File outputFilterSummaryFile = outputPrefix + ".filter_summary.json"
        File outputReportFile = outputPrefix + ".report.csv"
    }

    runtime {
        cpu: cores
        memory: memory
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        minPolyAlength: {description: "Minimum poly(A) tail length.", category: "advanced"}
        requirePolyA: {description: "Require FL reads to have a poly(A) tail and remove it.", category: "common"}
        logLevel: {description: "Set log level. Valid choices: (TRACE, DEBUG, INFO, WARN, FATAL).", category: "advanced"}
        inputBamFile: {description: "BAM input file.", category: "required"}
        primerFile: {description: "Barcode/primer fasta file.", category: "required"}
        outputPrefix: {description: "Output directory path + output file prefix.", category: "required"}
        cores: {description: "The number of cores to be used.", category: "advanced"}
        memory: {description: "The amount of memory available to the job.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        outputFLfile: {description: "Filtered reads output file."}
        outputFLindexFile: {description: "Index of filtered reads output file."}
        outputSTDERRfile: {description: "Refine STDERR log file."}
        outputConsensusReadsetFile: {description: "Refine consensus readset XML file."}
        outputFilterSummaryFile: {description: "Refine summary file."}
        outputReportFile: {description: "Refine report file."}
    }
}
