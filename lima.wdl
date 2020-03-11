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

task lima {
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
        ~{"--max-scored-barcode-pairs " + maxScoredBarcodePairs} \
        ~{"--max-scored-barcodes " + maxScoredBarcodes} \
        ~{"--max-scored-adapters " + maxScoredAdapters} \
        ~{"--min-passes " + minPasses} \
        ~{"--min-length " + minLength} \
        ~{"--max-input-length " + maxInputLength} \
        ~{"--min-ref-span " + minRefSpan} \
        ~{"--min-scoring-regions " + minScoringRegion} \
        ~{"--min-score " + minScore} \
        ~{"--min-end-score " + minEndScore} \
        ~{"--min-signal-increase " + minSignalIncrease} \
        ~{"--min-score-lead " + minScoreLead} \
        ~{true="--ccs" false="" ccsMode} \
        ~{true="--split-bam-named" false="" splitBamNamed} \
        ~{"--scored-adapter-ratio " + scoredAdapterRatio} \
        ~{"--peek " + peek} \
        ~{"--guess " + guess} \
        ~{"--guess-min-count " + guessMinCount} \
        ~{true="--peek-guess" false="" peekGuess} \
        ~{"--log-level " logLevel} \
        ~{"--num-threads " + cores} \
        ~{"--log-file " + outputPrefix + "_fl_stderr.log"} \
        ~{inputBamFile} \
        ~{barcodeFile} \
        ~{outputPrefix + ".fl.bam"}
    }

    output {
        File outputFLfile = outputPrefix + ".fl.bam"
        File outputFLindexFile = outputPrefix + ".fl.bam.pbi"
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
        libraryDesign: {description: "", category: ""}
        # outputs
    }
}
