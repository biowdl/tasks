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
        String outputDir
        String outputNamePrefix

        Int cores = 2
        String memory = "2G"
        Int timeMinutes = 180
        String dockerImage = "quay.io/biocontainers/isoseq3:3.3.0--0"
    }

    command {
        set -e
        mkdir -p "~{outputDir}"
        isoseq3 refine \
        --min-polya-length ~{minPolyAlength} \
        ~{true="--require-polya" false="" requirePolyA} \
        --log-level ~{logLevel} \
        --num-threads ~{cores} \
        --log-file "~{outputDir}/~{outputNamePrefix}.stderr.log" \
        ~{inputBamFile} \
        ~{primerFile} \
        "~{outputDir}/~{outputNamePrefix}.bam"
    }

    output {
        File outputFLNCfile = outputDir + "/" + outputNamePrefix + ".bam"
        File outputFLNCindexFile = outputDir + "/" + outputNamePrefix + ".bam.pbi"
        File outputConsensusReadsetFile = outputDir + "/" + outputNamePrefix + ".consensusreadset.xml"
        File outputFilterSummaryFile = outputDir + "/" + outputNamePrefix + ".filter_summary.json"
        File outputReportFile = outputDir + "/" + outputNamePrefix + ".report.csv"
        File outputSTDERRfile = outputDir + "/" + outputNamePrefix + ".stderr.log"
    }

    runtime {
        cpu: cores
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        minPolyAlength: {description: "Minimum poly(A) tail length.", category: "advanced"}
        requirePolyA: {description: "Require FL reads to have a poly(A) tail and remove it.", category: "common"}
        logLevel: {description: "Set log level. Valid choices: (TRACE, DEBUG, INFO, WARN, FATAL).", category: "advanced"}
        inputBamFile: {description: "BAM input file.", category: "required"}
        primerFile: {description: "Barcode/primer fasta file.", category: "required"}
        outputDir: {description: "Output directory path.", category: "required"}
        outputNamePrefix: {description: "Basename of the output files.", category: "required"}
        cores: {description: "The number of cores to be used.", category: "advanced"}
        memory: {description: "The amount of memory available to the job.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        outputFLNCfile: {description: "Filtered reads output file."}
        outputFLNCindexFile: {description: "Index of filtered reads output file."}
        outputSTDERRfile: {description: "Refine STDERR log file."}
        outputConsensusReadsetFile: {description: "Refine consensus readset XML file."}
        outputFilterSummaryFile: {description: "Refine summary file."}
        outputReportFile: {description: "Refine report file."}
    }
}
