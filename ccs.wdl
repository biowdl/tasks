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
        String logLevel = "WARN"
        Int minPasses = 3
        Int topPasses = 60
        Int minLength = 10
        Int maxLength = 50000
        Boolean byStrand = false
        Boolean skipPolish = false
        Boolean all = false
        Boolean subreadFallback = false
        Boolean allKinetics = false
        Boolean hifiKinetics = false
        Float minSnr = 2.5
        Float minReadQuality = 0.99

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
        --min-snr ~{minSnr} \
        --top-passes ~{topPasses} \
        --min-length ~{minLength} \
        --max-length ~{maxLength} \
        ~{true="--by-strand" false="" byStrand} \
        ~{true="--skip-polish" false="" skipPolish} \
        ~{true="--all" false="" all} \
        ~{true="--subread-fallback" false="" subreadFallback} \
        ~{true="--all-kinetics" false="" allKinetics} \
        ~{true="--hifi-kinetics" false="" hifiKinetics} \
        --min-rq ~{minReadQuality} \
        --log-level ~{logLevel} \
        --num-threads ~{threads} \
        ~{"--chunk " + chunkString} \
        ~{"--report-file " + outputPrefix + ".ccs_report.txt"} \
        ~{"--report-json " + outputPrefix + ".ccs.report.json"} \
        ~{"--log-file " + outputPrefix + ".ccs.stderr.log"} \
        ~{"--metrics-json " + outputPrefix + ".zmw_metrics.json.gz"} \
        ~{subreadsFile} \
        ~{outputPrefix + ".ccs.bam"}
    }

    output {
        File ccsBam = outputPrefix + ".ccs.bam"
        File ccsBamIndex = outputPrefix + ".ccs.bam.pbi"
        File ccsReport = outputPrefix + ".ccs_report.txt"
        File ccsJsonReport = outputPrefix + ".ccs.report.json"
        File ccsStderr = outputPrefix + ".ccs.stderr.log"
        File zmwMetrics = outputPrefix + ".zmw_metrics.json.gz"
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
        logLevel: {description: "Set log level. Valid choices: (TRACE, DEBUG, INFO, WARN, FATAL).", category: "advanced"}
        minPasses: {description: "Minimum number of full-length subreads required to generate ccs for a ZMW.", category: "advanced"}
        topPasses: {description: "Pick at maximum the top N passes for each ZMW.", category: "advanced"}
        minLength: {description: "Minimum draft length before polishing.", category: "advanced"}
        maxLength: {description: "Maximum draft length before polishing.", category: "advanced"}
        byStrand: {description: "Generate a consensus for each strand.", category: "advanced"}
        skipPolish: {description: "Only output the initial draft template (faster, less accurate).", category: "advanced"}
        all: {description: "Emit all ZMWs.", category: "advanced"}
        subreadFallback: {description: "Emit a representative subread, instead of the draft consensus, if polishing failed.", category: "advanced"}
        allKinetics: {description: "Calculate mean pulse widths (PW) and interpulse durations (IPD) for every ZMW.", category: "advanced"}
        hifiKinetics: {description: "Calculate mean pulse widths (PW) and interpulse durations (IPD) for every HiFi read.", category: "advanced"}
        minSnr: {description: "Minimum SNR of subreads to use for generating CCS.", category: "advanced"}
        minReadQuality: {description: "Minimum predicted accuracy in [0, 1].", category: "common"}
        subreadsIndexFile: {description: "Index for the subreads input file, required when using chunkString.", category: "advanced"}
        chunkString: {descpription: "Chunk string (e.g. 1/4, 5/5) for CCS.", category: "advanced"}
        threads: {description: "The number of threads to be used.", category: "advanced"}
        memory: {description: "The amount of memory available to the job.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        ccsBam: {description: "Consensus reads output file."}
        ccsBamIndex: {description: "Index of consensus reads output file."}
        ccsReport: {description: "Ccs report file."}
        ccsJsonReport: {description: "Ccs results json report file."}
        ccsStderr: {description: "Ccs STDERR log file."}
        zmwMetrics: {description: "ZMW metrics json file."}
    }
}
