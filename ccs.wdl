version 1.0

# Copyright (c) 2020 Leiden University Medical Center
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

task CCS {
    input {
        File subreadsFile
        String outputPrefix
        Int minPasses = 3
        Int minLength = 10
        Int maxLength = 50000
        Boolean byStrand = false
        Float minReadQuality = 0.99
        String logLevel = "WARN"

        File? subreadsIndexFile
        String? chunkString

        Int threads = 2
        String memory = "4G"
        Int timeMinutes = 1440
        String dockerImage = "quay.io/biocontainers/pbccs:6.0.0--h9ee0642_2"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPrefix})"
        ccs \
        --min-passes ~{minPasses} \
        --min-length ~{minLength} \
        --max-length ~{maxLength} \
        ~{true="--by-strand" false="" byStrand} \
        --min-rq ~{minReadQuality} \
        --log-level ~{logLevel} \
        --num-threads ~{threads} \
        ~{"--chunk " + chunkString} \
        ~{"--report-json " + outputPrefix + ".ccs.report.json"} \
        ~{"--log-file " + outputPrefix + ".ccs.stderr.log"} \
        ~{subreadsFile} \
        ~{outputPrefix + ".ccs.bam"}
    }

    output {
        File ccsBam = outputPrefix + ".ccs.bam"
        File ccsBamIndex = outputPrefix + ".ccs.bam.pbi"
        File ccsReport = outputPrefix + ".ccs.report.json"
        File ccsStderr = outputPrefix + ".ccs.stderr.log"
    }

    runtime {
        cpu: threads
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        subreadsFile: {description: "Subreads input file.", category: "required"}
        outputPrefix: {description: "Output directory path + output file prefix.", category: "required"}
        minPasses: {description: "Minimum number of full-length subreads required to generate ccs for a ZMW.", category: "advanced"}
        minLength: {description: "Minimum draft length before polishing.", category: "advanced"}
        maxLength: {description: "Maximum draft length before polishing.", category: "advanced"}
        byStrand: {description: "Generate a consensus for each strand.", category: "advanced"}
        minReadQuality: {description: "Minimum predicted accuracy in [0, 1].", category: "common"}
        logLevel: {description: "Set log level. Valid choices: (TRACE, DEBUG, INFO, WARN, FATAL).", category: "advanced"}
        subreadsIndexFile: {description: "Index for the subreads input file, required when using chunkString.", category: "advanced"}
        chunkString: {descpription: "Chunk string (e.g. 1/4, 5/5) for CCS.", category: "advanced"}
        threads: {description: "The number of threads to be used.", category: "advanced"}
        memory: {description: "The amount of memory available to the job.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        ccsBam: {description: "Consensus reads output file."}
        ccsBamIndex: {description: "Index of consensus reads output file."}
        ccsReport: {description: "Ccs results report file."}
        ccsStderr: {description: "Ccs STDERR log file."}
    }
}
