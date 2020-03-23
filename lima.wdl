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

task Lima {
    input {
        String libraryDesign = "same"
        Boolean scoreFullPass = false
        Int maxScoredBarcodePairs = 0
        Int maxScoredBarcodes = 0
        Int maxScoredAdapters = 0
        Int minPasses = 0
        Int minLength = 50
        Int maxInputLength = 0
        Float minRefSpan = 0.5
        Int minScoringRegion = 1
        Int minScore = 0
        Int minEndScore = 0
        Int minSignalIncrease = 10
        Int minScoreLead = 10
        Boolean ccsMode = false
        Boolean splitBamNamed = false
        Float scoredAdapterRatio = 0.25
        Int peek = 0
        Int guess = 0
        Int guessMinCount = 0
        Boolean peekGuess = false
        String logLevel = "WARN"
        File inputBamFile
        File barcodeFile
        String outputPrefix
        
        Int cores = 4
        String memory = "10G"
        String dockerImage = "quay.io/biocontainers/lima:1.11.0--0"
    }

    Map[String, String] libraryDesignOptions = {"same": "--same", "different": "--different", "neighbors": "--neighbors"}

    command {
        set -e
        mkdir -p "$(dirname ~{outputPrefix})"
        lima \
        ~{libraryDesignOptions[libraryDesign]} \
        ~{true="--score-full-pass" false="" scoreFullPass} \
        --max-scored-barcode-pairs ~{maxScoredBarcodePairs} \
        --max-scored-barcodes ~{maxScoredBarcodes} \
        --max-scored-adapters ~{maxScoredAdapters} \
        --min-passes ~{minPasses} \
        --min-length ~{minLength} \
        --max-input-length ~{maxInputLength} \
        --min-ref-span ~{minRefSpan} \
        --min-scoring-regions ~{minScoringRegion} \
        --min-score ~{minScore} \
        --min-end-score ~{minEndScore} \
        --min-signal-increase ~{minSignalIncrease} \
        --min-score-lead ~{minScoreLead} \
        ~{true="--ccs" false="" ccsMode} \
        ~{true="--split-bam-named" false="" splitBamNamed} \
        --scored-adapter-ratio ~{scoredAdapterRatio} \
        --peek ~{peek} \
        --guess ~{guess} \
        --guess-min-count ~{guessMinCount} \
        ~{true="--peek-guess" false="" peekGuess} \
        --log-level ~{logLevel} \
        --num-threads ~{cores} \
        ~{"--log-file " + outputPrefix + ".fl.stderr.log"} \
        ~{inputBamFile} \
        ~{barcodeFile} \
        ~{basename(outputPrefix) + ".fl.bam"}

        # Move commands below are needed because glob command does not find
        # multiple bam/bam.pbi/subreadset.xml files when not located in working
        # directory.
        mv "~{basename(outputPrefix)}.fl.json" "~{outputPrefix}.fl.json"
        mv "~{basename(outputPrefix)}.fl.lima.counts" "~{outputPrefix}.fl.lima.counts"
        mv "~{basename(outputPrefix)}.fl.lima.report" "~{outputPrefix}.fl.lima.report"
        mv "~{basename(outputPrefix)}.fl.lima.summary" "~{outputPrefix}.fl.lima.summary"
    }

    output {
        Array[File] outputFLfile = glob("~{basename(outputPrefix)}*.bam")
        Array[File] outputFLindexFile = glob("~{basename(outputPrefix)}*.bam.pbi")
        Array[File] outputFLxmlFile = glob("~{basename(outputPrefix)}*.subreadset.xml")
        File outputSTDERRfile = outputPrefix + ".fl.stderr.log"
        File outputJSONfile = outputPrefix + ".fl.json"
        File outputCountsFile = outputPrefix + ".fl.lima.counts"
        File outputReportFile = outputPrefix + ".fl.lima.report"
        File outputSummaryFile = outputPrefix + ".fl.lima.summary"
    }

    runtime {
        cpu: cores
        memory: memory
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        libraryDesign: {description: "Barcode structure of the library design.", category: "required"}
        scoreFullPass: {description: "Only use subreads flanked by adapters for barcode identification.", category: "advanced"}
        maxScoredBarcodePairs: {description: "Only use up to N barcode pair regions to find the barcode, 0 means use all.", category: "advanced"}
        maxScoredBarcodes: {description: "Analyze at maximum the provided number of barcodes per ZMW, 0 means deactivated.", category: "advanced"}
        maxScoredAdapters: {description: "Analyze at maximum the provided number of adapters per ZMW, 0 means deactivated.", category: "advanced"}
        minPasses: {description: "Minimal number of full passes.", category: "common"}
        minLength: {description: "Minimum sequence length after clipping.", category: "common"}
        maxInputLength: {description: "Maximum input sequence length, 0 means deactivated.", category: "advanced"}
        minRefSpan: {description: "Minimum reference span relative to the barcode length.", category: "advanced"}
        minScoringRegion: {description: "Minimum number of barcode regions with sufficient relative span to the barcode length.", category: "advanced"}
        minScore: {description: "Reads below the minimum barcode score are removed from downstream analysis.", category: "common"}
        minEndScore: {description: "Minimum end barcode score threshold is applied to the individual leading and trailing ends.", category: "advanced"}
        minSignalIncrease: {description: "The minimal score difference, between first and combined, required to call a barcode pair different.", category: "advanced"}
        minScoreLead: {description: "The minimal score lead required to call a barcode pair significant.", category: "common"}
        ccsMode: {description: "CCS mode, use optimal alignment options.", category: "common"}
        splitBamNamed: {description: "Split BAM output by resolved barcode pair name.", category: "common"}
        scoredAdapterRatio: {description: "Minimum ratio of scored vs sequenced adapters.", category: "advanced"}
        peek: {description: "Demux the first N ZMWs and return the mean score, 0 means peeking deactivated.", category: "advanced"}
        guess: {description: "Try to guess the used barcodes, using the provided mean score threshold, 0 means guessing deactivated.", category: "advanced"}
        guessMinCount: {description: "Minimum number of ZMWs observed to whitelist barcodes.", category: "advanced"}
        peekGuess: {description: "Try to infer the used barcodes subset, by peeking at the first 50,000 ZMWs.", category: "advanced"}
        logLevel: {description: "Set log level. Valid choices: (TRACE, DEBUG, INFO, WARN, FATAL).", category: "advanced"}
        inputBamFile: {description: "BAM input file.", category: "required"}
        barcodeFile: {description: "Barcode/primer fasta file.", category: "required"}
        outputPrefix: {description: "Output directory path + output file prefix.", category: "required"}
        cores: {description: "The number of cores to be used.", category: "advanced"}
        memory: {description: "The amount of memory available to the job.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        outputFLfile: {description: "Demultiplexed reads output file."}
        outputFLindexFile: {description: "Index of demultiplexed reads output file."}
        outputSTDERRfile: {description: "Lima STDERR log file."}
        outputJSONfile: {description: "Lima JSON file."}
        outputCountsFile: {description: "Lima counts file."}
        outputReportFile: {description: "Lima report file."}
        outputSummaryFile: {description: "Lima summary file."}
    }
}
