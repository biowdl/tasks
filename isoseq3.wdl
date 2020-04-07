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
        Int timeMinutes = ceil(size(inputBamFile, "G") * 240 / cores)
    }

    command <<<
        set -e
        mkdir -p "$(dirname ~{outputPrefix})"

        # Create a unique output name base on the input bam file.
        bamBasename="$(basename ~{inputBamFile})"
        bamNewName="${bamBasename/fl/flnc}"
        folderDirname="$(dirname ~{outputPrefix})"
        combinedOutput="${folderDirname}/${bamNewName}"

        isoseq3 refine \
        --min-polya-length ~{minPolyAlength} \
        ~{true="--require-polya" false="" requirePolyA} \
        --log-level ~{logLevel} \
        --num-threads ~{cores} \
        --log-file "${bamNewName}.stderr.log" \
        ~{inputBamFile} \
        ~{primerFile} \
        ${bamNewName}

        # Copy commands below are needed because naming schema for Refine output
        # can not be correctly handled in the WDL output section.
        cp "${bamNewName}" "${combinedOutput}"
        cp "${bamNewName}.pbi" "${combinedOutput}.pbi"
        cp "${bamNewName/bam/consensusreadset}.xml" "${combinedOutput/bam/consensusreadset}.xml"
        cp "${bamNewName/bam/filter_summary}.json" "${combinedOutput/bam/filter_summary}.json"
        cp "${bamNewName/bam/report}.csv" "${combinedOutput/bam/report}.csv"
        cp "${bamNewName}.stderr.log" "${combinedOutput}.stderr.log"
    >>>

    output {
        Array[File] outputFLNCfile = glob("*.bam")
        Array[File] outputFLNCindexFile = glob("*.bam.pbi")
        Array[File] outputConsensusReadsetFile = glob("*.consensusreadset.xml")
        Array[File] outputFilterSummaryFile = glob("*.filter_summary.json")
        Array[File] outputReportFile = glob("*.report.csv")
        Array[File] outputSTDERRfile = glob("*.stderr.log")
    }

    runtime {
        cpu: cores
        memory: memory
        docker: dockerImage
        time_minutes: timeMinutes
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
        timeMinutes: {description: "The time (in minutes) it will take for this task to complete.", category: "advanced"}

        # outputs
        outputFLNCfile: {description: "Filtered reads output file."}
        outputFLNCindexFile: {description: "Index of filtered reads output file."}
        outputSTDERRfile: {description: "Refine STDERR log file."}
        outputConsensusReadsetFile: {description: "Refine consensus readset XML file."}
        outputFilterSummaryFile: {description: "Refine summary file."}
        outputReportFile: {description: "Refine report file."}
    }
}
