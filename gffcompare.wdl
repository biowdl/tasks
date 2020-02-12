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

task GffCompare {
    input {
        File? inputGtfList
        Array[File] inputGtfFiles
        File referenceAnnotation
        String? outputDir
        String outPrefix = "gffcmp" # gffcmp is the default used by the program as well. This
        # needs to be defined in order for the output values to be consistent and correct.
        File? genomeSequences
        Int? maxDistanceFreeEndsTerminalExons
        Int? maxDistanceGroupingTranscriptStartSites
        String? namePrefix
        Boolean C = false
        Boolean A = false
        Boolean X = false
        Boolean K = false
        Boolean snCorrection = false
        Boolean precisionCorrection = false
        Boolean discardSingleExonTransfragsAndReferenceTranscripts = false
        Boolean discardSingleExonReferenceTranscripts = false
        Boolean noTmap = false
        Boolean verbose = false
        Boolean debugMode = false

        String dockerImage = "quay.io/biocontainers/gffcompare:0.10.6--h2d50403_0"

        # This workaround only works in the input section.
        # Issue addressed at https://github.com/openwdl/wdl/pull/263
        File? noneFile # This is a wdl workaround. Please do not assign!
    }

    # This allows for the creation of output directories
    String dirPrefix = if defined(outputDir)
        then select_first([outputDir]) + "/"
        else ""
    String totalPrefix = dirPrefix + outPrefix

    command {
        set -e
        ~{"mkdir -p " + outputDir}
        gffcompare \
        -r ~{referenceAnnotation} \
        ~{"-o '" + totalPrefix + "'"} \
        ~{"-s " + genomeSequences} \
        ~{"-e " + maxDistanceFreeEndsTerminalExons} \
        ~{"-d " + maxDistanceGroupingTranscriptStartSites} \
        ~{"-p " + namePrefix} \
        ~{true="-C" false="" C} \
        ~{true="-A" false="" A} \
        ~{true="-X" false="" X} \
        ~{true="-K" false="" K} \
        ~{true="-R" false="" snCorrection} \
        ~{true="-Q" false="" precisionCorrection} \
        ~{true="-M" false="" discardSingleExonTransfragsAndReferenceTranscripts} \
        ~{true="-N" false="" discardSingleExonReferenceTranscripts} \
        ~{true="-T" false="" noTmap} \
        ~{true="-V" false="" verbose} \
        ~{true="D" false="" debugMode} \
        ~{"-i " + inputGtfList} \
        ~{sep=" " inputGtfFiles}
    }

    # Output of gffcompare is not stable. It depends on the number of files in the input.
    Int noFilesGtfList = if defined(inputGtfList)
        then length(read_lines(select_first([inputGtfList])))
        else 0
    Int noInputFiles = length(inputGtfFiles)
    Boolean oneFile = (noFilesGtfList + noInputFiles) == 1
    String annotatedName = if oneFile
        then "annotated"
        else "combined"

    # Check if a redundant .gtf will be created
    Boolean createRedundant = C || A || X

    output {
        File annotated = totalPrefix + "." + annotatedName + ".gtf"
        File loci = totalPrefix + ".loci"
        File stats = totalPrefix + ".stats"
        File tracking = totalPrefix + ".tracking"
        # noneFile is not stable. Please replace this as soon as wdl spec allows
        File? redundant = if createRedundant
            then totalPrefix + ".redundant.gtf"
            else noneFile
        File? missedIntrons = if debugMode
            then totalPrefix + ".missed_introns.gtf"
            else noneFile
    }

    runtime {
       docker: dockerImage
    }

    parameter_meta {
        inputGtfList: {description: "Equivalent to gffcompare's `-i` option.", category: "advanced"}
        inputGtfFiles: {description: "The input GTF files.", category: "required"}
        referenceAnnotation: {description: "The GTF file to compare with.", category: "required"}
        outputDir: {description: "The location the output should be written.", category: "common"}
        outPrefix: {description: "The prefix for the output.", category: "advanced"}
        genomeSequences: {description: "Equivalent to gffcompare's `-s` option.", category: "advanced"}
        maxDistanceFreeEndsTerminalExons: {description: "Equivalent to gffcompare's `-e` option.", category: "advanced"}
        maxDistanceGroupingTranscriptStartSites: {description: "Equivalent to gffcompare's `-d` option.", category: "advanced"}
        namePrefix: {description: "Equivalent to gffcompare's `-p` option.", category: "advanced"}
        C: {description: "Equivalent to gffcompare's `-C` flag.", category: "advanced"}
        A: {description: "Equivalent to gffcompare's `-A` flag.", category: "advanced"}
        X: {description: "Equivalent to gffcompare's `-X` flag.", category: "advanced"}
        K: {description: "Equivalent to gffcompare's `-K` flag.", category: "advanced"}
        snCorrection: {description: "Equivalent to gffcompare's `-R` flag.", category: "advanced"}
        precisionCorrection: {description: "Equivalent to gffcompare's `-Q` flag.", category: "advanced"}
        discardSingleExonTransfragsAndReferenceTranscripts: {description: "Equivalent to gffcompare's `-M` flag.", category: "advanced"}
        discardSingleExonReferenceTranscripts: {description: "Equivalent to gffcompare's `-N` flag.", category: "advanced"}
        noTmap: {description: "Equivalent to gffcompare's `-T` flag.", category: "advanced"}
        verbose: {description: "Equivalent to gffcompare's `-V` flag.", category: "advanced"}
        debugMode: {description: "Equivalent to gffcompare's `-D` flag.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }

    meta {
        WDL_AID: {
            exclude: ["noneFile"]
        }
    }
}